import argparse

parser = argparse.ArgumentParser(
    prog = 'Make Dataset',
    description = 'This cleanes a bunch of pdb files or just one.',
)
parser.add_argument('-i', '--input_dir')
parser.add_argument('-o', '--output_dir')

import os
import pandas as pd
from biopandas.pdb import PandasPdb
from alive_progress import alive_bar
import numpy as np

def save_to_pdb(structure:pd.DataFrame, path:str)->None:
    """
    Save extracted segments to a pdb file
    """
    pdb_saver = PandasPdb()
    pdb_saver.df["ATOM"] = structure
    pdb_saver.to_pdb(path, records = ["ATOM"])

def clean_pdb(file:str, output_file:str):
    if file.endswith(".pdb"):
        pdb_df = PandasPdb().read_pdb(file)
        pdb_df.label_models()
        pdb_df = pdb_df.get_model(1)
        pdb = pdb_df.df["ATOM"]
        pdb.drop(columns=["model_id"], inplace=True)
        mask = pdb["alt_loc"] == ""
        pdb = pdb[mask].copy()
        mask = pdb["element_symbol"] != "H"
        pdb = pdb[mask].copy()
        mask = pdb["insertion"] == ""
        pdb = pdb[mask].copy()

        chains = pdb["chain_id"].unique()

        pdb["nresidue_number"] = np.nan

        for chain in chains:
            resnums_to_new_resnums = {old_rn:i+1 for i, old_rn in enumerate(pdb.loc[pdb["chain_id"] == chain]["residue_number"].unique()) if old_rn != i+1}
            for (o, n) in resnums_to_new_resnums.items():
                pdb.loc[(pdb["chain_id"] == chain)  & (pdb["residue_number"]==o), "nresidue_number"] = n
        pdb["residue_number"] = pdb.apply(lambda x: x["nresidue_number"] if not np.isnan(x["nresidue_number"]) else x["residue_number"], axis=1)
        pdb["residue_number"]= pdb["residue_number"].astype("int32")
        pdb.drop(columns=["nresidue_number"], inplace=True)
        pdb.reset_index(inplace=True, drop=True)
        pdb["atom_number"] = pdb.index+1
        save_to_pdb(pdb, os.path.join(output_file, file))



def clean_pdbs(input_dir:str, output_dir:str):
    with alive_bar(len(os.listdir(input_dir)), force_tty=True) as bar:
        for file in os.listdir(input_dir):
            if file.endswith(".pdb"):
                pdb_df = PandasPdb().read_pdb(os.path.join(input_dir, file))
                pdb_df.label_models()
                pdb_df = pdb_df.get_model(1)
                pdb = pdb_df.df["ATOM"]
                pdb.drop(columns=["model_id"], inplace=True)
                mask = pdb["alt_loc"] == ""
                pdb = pdb[mask].copy()
                mask = pdb["element_symbol"] != "H"
                pdb = pdb[mask].copy()
                mask = pdb["insertion"] == ""
                pdb = pdb[mask].copy()

                chains = pdb["chain_id"].unique()

                pdb["nresidue_number"] = np.nan

                for chain in chains:
                    resnums_to_new_resnums = {old_rn:i+1 for i, old_rn in enumerate(pdb.loc[pdb["chain_id"] == chain]["residue_number"].unique()) if old_rn != i+1}
                    for (o, n) in resnums_to_new_resnums.items():
                        pdb.loc[(pdb["chain_id"] == chain)  & (pdb["residue_number"]==o), "nresidue_number"] = n
                pdb["residue_number"] = pdb.apply(lambda x: x["nresidue_number"] if not np.isnan(x["nresidue_number"]) else x["residue_number"], axis=1)
                pdb["residue_number"]= pdb["residue_number"].astype("int32")
                pdb.drop(columns=["nresidue_number"], inplace=True)
                pdb.reset_index(inplace=True, drop=True)
                pdb["atom_number"] = pdb.index+1
                save_to_pdb(pdb, os.path.join(output_dir, file))
            bar()

if __name__ == "__main__":

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    if input_dir == output_dir:
        raise ValueError("The input folder and the output folder shouldn't be the same one due to safety reasons.")
    if os.path.isdir(input_dir) and os.path.isdir(output_dir):
        clean_pdbs(input_dir, output_dir)
    elif os.path.isfile(input_dir) and os.path.isfile(output_dir):
        clean_pdb(input_dir, output_dir)
    elif os.path.isfile(input_dir) and os.path.isdir(output_dir):
        head, tail = os.path.split(input_dir)
        clean_pdb(input_dir, os.path.join(output_dir, tail))
    else:
        print("You did something wrong. Try again.")

