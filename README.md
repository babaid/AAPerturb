# AAPerturb
A C++ library for the creation of a large dataset of amino acid sidechain perturbations, own PDB Parser code included and some other things related as partly described in [Pre-training of Graph Neural Network for Modeling Effects of Mutations on Protein-Protein Binding Affinity](https://arxiv.org/abs/2008.12473).

This is a faster implementation of my AA-Perturbation library in C++ instead of python.

## PDB Parser
It has an own PDB parser which is mostly educational/for myself as there is probably something available that does this.

The proteins are basically represented as a map chainID -> Vector of Residues, where the Residues contain the atoms.
A particular reason I implemented this myself, is that most libraries that are available in python overcomplexify the task of parsing PDBs,
introducing weird classes with crazy functionality that is useful in many cases but in my case I just want to do geometric transormations on simple ATOM records.

For some of the implementations refer to Biopython, rdkit, and the really useful Biopandas.


## Geometrical operations

I would say my code is anything else than highly optimized. There are different ways to implement rotations, translations and calculations of distances for molecules.
Although it would have been more efficient to work with matrices of atomic coordinates, this would mean to keep track of exactly how those matrices are ordered in terms of atoms and residues.
So to keep it simple, and bookkeeping reasons I decided to work with single coordinates of the atoms. Each split into residues and chains. If someone has a more efficient idea please contact me.

## Random Perturbations

The main goal of the package is to perform random perturbations the sidechain of a random amino acid, this amino acid should reside
on the interface between two chains of the protein-protein complex.
The created data set can be used after for an autoencoder-like machine learning approach to capture PPI's and effects of mutations proteins.

The execution flow is as follows:
1. Find interface residues on the PP complex, given a cutoff value
2. Choose random interface residue
3. Perturb chosen residue

For now this perturbation is just torsion about the sidechain axes which could in theory freely rotate.
In my head the only condition for acceptance of a conformation after a perturbation, is that there are no clashes between atoms.
It is also possible to sample the perturbation \Chi angles from physically relevant distributions as defined for example in [Dunbrack (2011)](http://dunbrack.fccc.edu/lab/bbdep2010).


Finally, to be specific about what I would like to implement in the near future, as some extra functionality.
Just like in these geometrical transformations are perturbations, we can look at mutations in the same way. While the concept proposed in the first mentioned paper is to recreate atomic coordinates after a perturbation, a similar operation could be used in an alchemical way, recreating both coordinates, and adding atoms/graph nodes to a molecular graph, with the goal to reduce its strain/energy or recreate the original structure.
In principle this could be used to predict the probability of certain mutations, which is obviously something cool and useful.


## PP Interfaces

One more thing that I havent found yet, is a program that finds the interfaces between PP complexes, so I provide it right away.


## Requirements

The library was built using C++23, with g++ and gcc version 13.
To build it clone this project and inside the project directory:


```
mkdir build
cd build
cmake ..
make -jN aaperturb
```

or if you want to use the interface residue finder executable:

```
make -jN iffind
```

### Cleaning of PDB files.

To produce output that is useful, a preprocessing step may be required, which deals with removal of waters, deprotonation, removal of alternate locations, removal of insertions and reindexing atoms and residues.
This is crucial to the main program to run without bugs and random errors, and also the less atoms there are the less expensive calculations get.
The residue numbering has to start in each chain at 1.
The part of cleaning the PDBs seemed easier to do with python with the already available PandasPDB package. The script is located at scripts/pdbcleaner.py.
It needs following packages: biopandas, pandas, numpy, alive_progress
You can run it as follows:
```
python pdbcleaner.py -i [INPUT_DIR/FILE] -o [OUTPUT_DIR/FILE]
```

It is going to clean all the files and save them. Specifying the input dir to be the same as the output dir is prohibited due to the safety of your dataset.








