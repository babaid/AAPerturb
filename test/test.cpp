//
// Created by babaid on 11.09.23.
//
#include<map>
#include<vector>
#include<argparse/argparse.hpp>
#include<filesystem>
#include<format>
#include<iostream>
#include "molecules.h"
namespace fs = std::filesystem;
bool verbose=false;

int main(int argc, char *argv[]) {

    argparse::ArgumentParser program("aaperturb-singlerot");
    program.add_argument("-v", "--verbose")
            .help("Enable verbose mode. It is useful if you want to know if something is happening.")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("-i", "--input-file")
            .required()
            .help("The input PDB file which we will change.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output files.");
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }


    if(program["-v"] == true)
    {
        std::cout << "Running in verbose mode." << std::endl;
        verbose = true;
    }
    fs::path input_file;
    if (auto ifn = program.present("-i"))
    {
        input_file = *ifn;
        if(fs::is_directory(input_file))
        {
            std::cout << "Input should be a valid path to a PDB file. Exiting..." << std::endl;
            std::exit(1);
        }
    }
    fs::path output_dir;
    if (auto ofn = program.present("-o"))
    {
        output_dir = *ofn;
        if(!fs::is_directory(output_dir))
        {
            std::cout << "The directory provided does not exist. Creating directory " << output_dir << std::endl;
            fs::create_directory(output_dir);
        }
    }

    RandomPerturbator structure(input_file, verbose);

    structure.calculateDistanceMatrix();
    structure.saveDistMat(output_dir/"distma.tsv");
    std::cout << "You can choose from the following interface resiudes: " << std::endl;


    char chain{'0'};
    int resid{-1};
    Residue refres;
    while ( (chain == '0' && resid == -1 ) && refres.atom_coords.empty())
    {

        std::cout << "Choose one of the residues!\nChain ID: ";
        std::cin >> chain;
        std::cout << "Residue Index: ";
        std::cin >> resid;
        refres = structure.getResidue(chain, resid);
        if (refres.resName == "GLY" || refres.resName == "PRO" || refres.resName == "ALA") {
            std::cout << "Chosen residue is: " << refres.resName << ", which cant be rotated around any axis."
                      << std::endl;
            chain = '0';
            resid = -1;
            continue;
        }

    }

    auto rmsd = structure.rotateResidueSideChain(chain, resid);
    rmsd = structure.rotateResidueAroundBackbone(chain, resid);
    std::vector<std::string> comments;
    std::string comment1 = std::format("PERTURBATION: /{}:{}", chain, std::to_string(resid));
    std::string comment2 = std::format("RMSD: {}", rmsd);
    comments.push_back(comment1);
    comments.push_back(comment2);
    structure.savePDBStructure(output_dir/"rotated.pdb", comments);


}


