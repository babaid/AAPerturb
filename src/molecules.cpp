//
// Created by babaid on 09.09.23.
//
#include<iostream>
#include<algorithm>
#include<string>
#include<fstream>
#include<filesystem>
#include "molecules.h"
#include "pdbparser.h"

namespace fs=std::filesystem;


bool Atom::operator==(const Atom &atom) const {
    bool ser = this->serial == atom.serial;
    bool name = this->name == atom.name;
    bool rname = this->resName == atom.resName;
    bool el = this->element == atom.element;
    return ser && name && rname && el;
}



PDBStructure::PDBStructure(const fs::path path, bool exclude_waters, bool, bool){


    chains = std::make_unique<std::map<char, std::vector<Residue>>>();
    numResiduesPerChain = std::make_unique<std::map<char, unsigned int>>();
    std::ifstream pdbFile(path);

    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << path << std::endl;
    }

    std::string line;

    unsigned long long numAtoms{0}; 
    while (std::getline(pdbFile, line)) {
        if (line.compare(0, 4, "ATOM") == 0 || ((line.compare(0, 6, "HETATM") == 0) && !exclude_waters) )  {
            Atom atom;
            atom.serial = std::stoi(line.substr(6, 5));
            atom.name = line.substr(12, 4);
            atom.altLoc = line[16];
            atom.resName = line.substr(17, 3);
            atom.chainID = line[21];
            atom.resSeq = std::stoi(line.substr(22, 4))-1;
            atom.coords= {std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)), std::stod(line.substr(46, 8))};
            atom.occupancy = std::stod(line.substr(54, 6));
            atom.tempFactor = std::stod(line.substr(60, 6));
            atom.element = line.substr(76, 2);

            // Remove whitespace from the atom name
            atom.name.erase(std::remove_if(atom.name.begin(), atom.name.end(), ::isspace), atom.name.end());

            // Check if this chain is already in the map
            if (chains->find(atom.chainID) == chains->end()) {
                (*chains)[atom.chainID] = std::vector<Residue>();
                (*numResiduesPerChain)[atom.chainID] = 0;
            }

            // Check if this residue is already in the chain's residues
            bool found = false;
            for (auto& residue : chains->at(atom.chainID)) {
                if (residue.resSeq == atom.resSeq && residue.resName == atom.resName) {
                    residue.atoms.emplace_back(std::move(atom));
                    found = true;
                    break;
                }
            }

            // If the residue doesn't exist, create a new one
            if (!found) {
                Residue newResidue;
                newResidue.chainID = atom.chainID;
                newResidue.resSeq = atom.resSeq; //residueCounter;
                newResidue.resName = atom.resName;
                newResidue.atoms.emplace_back(std::move(atom));
                chains->at(atom.chainID).emplace_back(std::move(newResidue));
                (*numResiduesPerChain)[atom.chainID]++;
             
            }
            numAtoms++;
        }
    }

    pdbFile.close();

}
