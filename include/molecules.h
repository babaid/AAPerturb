//Header for molecule data structures
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB__MOLECULES_H
#define AAPERTURB__MOLECULES_H
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

//A class for later, not used RN but certainly cool
class PDBStructure{
    private:
        bool verbose{true};
        std::map<char, std::vector<Residue>> chains; //The actual structure
        std::map<char, std::vector<unsigned>> interfaceResidues;

        unsigned int numAtoms{0}; // Number of atoms.

        std::vector<std::vector<double>> dist_mat;//(numAtoms, std::vector<double>(numAtoms, 0.0));
        std::vector<std::vector<double>> atomic_feat_mat; //Only has to be calculated once.
    public:
        PDBStructure(const fs::path, bool);

        void calculateDistanceMatrix(); //Calculates distance matrix of one protein
        void findInterfaceResidues(double);


        std::vector<std::vector<double>> calculateLocalDistanceMatrix(Residue& refres); //you can use any reference residue. This is extremely useful
        void updateDistanceMatrixLocally(std::vector<std::vector<double>>& newPart, char chain, unsigned int resNum); //subsitutes part of the distance matrix. Returns a new one. We dont want to change the original one.

        [[nodiscard]] std::pair<char , std::size_t> chooseRandomResidue() const;
        double rotateResidueSidechainRandomly(char, std::size_t);

        double rotateResidueSideChain(char, std::size_t );

        double rotateResidueAroundBackbone(char, std::size_t);

        void calculateAtomicFeatureMatrix();
        void savePDBStructure(fs::path, const std::vector<std::string>&);
        void saveCoords(fs::path);
        void saveFeatureMat(fs::path, bool);
        void saveDistMat(fs::path);
        void saveLocalDistMat(fs::path, char, int);
        void getNumberOfResidues();
        void getInterfaceResidues();
        std::vector<std::vector<double>> getDistMat();
        Residue getResidue(char, unsigned);
        void setResidue(const Residue&);

    };

class Protein: public PDBStructure
{

};
//std::vector<std::vector<double>> updateDistanceMatrixLocally(const PDBStructure& ,std::vector<std::vector<double>>& newPart, char chain, unsigned int resNum);
std::vector<double> one_hot(unsigned index, unsigned num_classes);

#endif //AA_SIDECHAIN_PERTURBATION_CPP_MOLECULES_H


