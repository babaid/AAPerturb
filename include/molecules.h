//Header for molecule data structures
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB_MOLECULES_H
#define AAPERTURB_MOLECULES_H
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
    bool operator==(const Atom& at) const;
};

double calculateDistance(const Atom& atom1, const Atom& atom2);

//A whole residue of atoms in a PDB file
struct Residue {
    char chainID;
    int resSeq;
    std::string resName;
    std::vector<Atom> atoms;
    std::vector<std::string> atom_seq;
    std::vector<std::vector<double>> atom_coords;
};


struct Protein
{
    std::vector<char> chain_ordering;
    std::map<char, std::vector<Residue>> chains;
    unsigned int numAtoms{0};
    explicit Protein(fs::path&);
    void saveToPDB(fs::path&, const std::vector<std::string>&);

};

/*
 * This class takes a Protein in PDB file format as input and allows to create other classes that manipulate the protein structure.
 */
class RandomPerturbator
{
    bool verbose;
    Protein protein;
    std::map<char, std::vector<unsigned>> interfaceResidues;
    double maxRotAngleBB, maxRotAngleSCH;

    //std::vector<std::vector<double>> dist_mat;

public:
    RandomPerturbator(fs::path&, bool);

    //void calculateDistanceMatrix(); //Calculates distance matrix of one protein
    std::vector<std::vector<double>> calculateLocalDistanceMatrix(Residue& refres); //you can use any reference residue. This is extremely useful
    double calculateRMSD(const Residue& ref_res);

    //also magic
    void findInterfaceResidues(double);
    void saveInterfaceResidues(fs::path&);
    [[nodiscard]] std::pair<char , std::size_t> chooseRandomResidue() const;



    //this is where the magic is happening
    void rotateResidueSidechainRandomly(char, std::size_t );
    void rotateResidueAroundBackboneRandomly(char, std::size_t);


    //savers
    //void saveCoords(const fs::path&);
    //void saveDistMat(const fs::path&);
    //void saveLocalDistMat(const fs::path&, char, int);
    void saveToPDB(fs::path&, const std::vector<std::string>&);

    //getters
    void getNumberOfResiduesPerChain();
    void getInterfaceResidues();
    std::vector<std::vector<double>> getDistMat();
    Residue getResidue(char, unsigned);



    //setters
    void setResidue(const Residue&);
    void setMaxRotAngleBB(double);
    void setMaxRotAngleSCH(double);
    //void setDistanceMatrixLocally(std::vector<std::vector<double>>& newPart, char chain, unsigned int resNum); //subsitutes part of the distance matrix. Returns a new one. We dont want to change the original one.

};


#endif //AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H


