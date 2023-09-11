//Header for molecule data structures
// Created by babaid on 09.09.23.
//

#ifndef AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H
#define AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H
#include<vector>
#include<string>
#include<valarray>
#include<map>
#include<filesystem>
namespace fs = std::filesystem;
//An atom record from a PDB file
struct Atom {
    int serial;
    int resSeq;
    char altLoc;
    char chainID;
    double occupancy;
    double tempFactor;
    std::valarray<double> coords{0., 0., 0};
    std::string name;
    std::string resName;
    std::string element;

};

//A whole residue of atoms in a PDB file
struct Residue {
    char chainID;
    int resSeq;
    std::string resName;
    std::vector<Atom> atoms;
    std::vector<std::string> atom_seq;
    std::vector<std::vector<double>> atom_coords;

};


//A class for later, not used RN but certainly cool
/*class Protein{
    private:
        std::map<char, std::vector<Residue>> chains;
        unsigned long numAtoms;

    public:
        Protein(fs::path);
        Protein(std::string);

    };*/

#endif //AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H
