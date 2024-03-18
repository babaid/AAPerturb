//Header for molecule data structures
#ifndef AAPERTURB_MOLECULES_H
#define AAPERTURB_MOLECULES_H
#include "spdlog/spdlog.h"
#include<vector>
#include<string>
#include<array>
#include<map>
#include<filesystem>

namespace fs = std::filesystem;
using Vector3 = std::array<double, 3>;

//An atom record from a PDB file
struct Atom {
    size_t serial;
    size_t resSeq;
    char altLoc;
    char chainID;
    double occupancy;
    double tempFactor;
    Vector3 coords{0., 0., 0};
    std::string name;
    std::string resName;
    std::string element;
    bool operator==(const Atom& at) const;
};

double calculateDistance(const Atom& atom1, const Atom& atom2);

//A whole residue of atoms in a PDB file
struct Residue {
    char chainID;
    size_t resSeq;
    std::string resName;
    std::vector<Atom> atoms;
    std::vector<std::string> atom_seq;
    std::vector<std::vector<double>> atom_coords;
};


struct Protein
{
    std::vector<char> chain_ordering;
    std::map<char, std::vector<Residue>> chains;
    size_t numAtoms{0};
    explicit Protein(fs::path&);
    void saveToPDB(fs::path&, const std::vector<std::string>&);

};

/*
 * This class takes a Protein in PDB file format as input and allows to create other classes that manipulate the protein structure.
 */
class RandomPerturbator
{
    std::shared_ptr<spdlog::logger> logger;
    Protein protein;
    std::map<char, std::vector<size_t>> interfaceResidues;
    double maxRotAngleBB, maxRotAngleSCH;

public:
    explicit RandomPerturbator(fs::path&);
    explicit RandomPerturbator(fs::path&, double maxRotationBB, double maxRotationSCH);

    std::vector<std::vector<double>> calculateLocalDistanceMatrix(Residue& refres); //you can use any reference residue. This is extremely useful
    double calculateRMSD(const Residue& ref_res);

    void findInterfaceResidues(double);
    void saveInterfaceResidues(fs::path&);

    [[nodiscard]] std::pair<char , std::size_t> chooseRandomResidue() const;


    //this is where the magic is happening
    void rotateResidueSidechainRandomly(char, std::size_t );
    void rotateResidueAroundBackboneRandomly(char, std::size_t);

    void saveToPDB(fs::path&, const std::vector<std::string>&);


    //getters
    void getNumberOfResiduesPerChain();
    std::map<char, std::vector<size_t>> getInterfaceResidues();

    Residue getResidue(char, size_t);
    void printInterfaceResidues();

    //setters
    void setResidue(const Residue&);
    void setMaxRotAngleBB(double);
    void setMaxRotAngleSCH(double);
};


#endif //AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H


