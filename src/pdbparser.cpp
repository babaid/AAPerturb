#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include<algorithm>
#include<random>
#include<iomanip>

#include "../include/pdbparser.h"
#include "../include/geometry.h"
#include "../include/constants.h"

namespace fs = std::filesystem;

std::map<char, std::vector<Residue*>> parsePDB(const  fs::path& filename, bool excludewaters, bool deprotonate) {
    std::map<char, std::vector<Residue*>>* chainMap = new std::map<char, std::vector<Residue*>>();
    std::ifstream pdbFile(filename);

    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return *chainMap;
    }

    std::string line;

    while (std::getline(pdbFile, line)) {
        if (line.compare(0, 4, "ATOM") == 0 || ((line.compare(0, 6, "HETATM") == 0) && !excludewaters) )  {
            Atom* atom = new Atom();
            atom->serial = std::stoi(line.substr(6, 5));
            atom->name = line.substr(12, 4);
            atom->altLoc = line[16];
            atom->resName = line.substr(17, 3);
            atom->chainID = line[21];
            atom->resSeq = std::stoi(line.substr(22, 4));
            atom->coords= {std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)), std::stod(line.substr(46, 8))};
            atom->occupancy = std::stod(line.substr(54, 6));
            atom->tempFactor = std::stod(line.substr(60, 6));
            atom->element = line.substr(76, 2);

            // Remove whitespace from the atom name
            atom->name.erase(std::remove_if(atom->name.begin(), atom->name.end(), ::isspace), atom->name.end());

            // Check if this chain is already in the map
            if (chainMap->find(atom->chainID) == chainMap->end()) {
                (*chainMap)[atom->chainID] = std::vector<Residue*>();
            }

            // Check if this residue is already in the chain's residues
            bool found = false;
            for (Residue* residue : (*chainMap)[atom->chainID]) {
                if (residue->resSeq == atom->resSeq && residue->resName == atom->resName) {
                    if(atom->element == "H" && deprotonate) continue;
                    else residue->atoms.push_back(*atom);
                    found = true;
                    break;
                }
            }

            // If the residue doesn't exist, create a new one
            if (!found) {
                Residue* newResidue = new Residue();
                newResidue->chainID = atom->chainID;
                newResidue->resSeq = atom->resSeq;
                newResidue->resName = atom->resName;
                if(atom->element == "H" && deprotonate) continue;
                else newResidue->atoms.push_back(*atom);
                //newResidue.atom_coords.push_back({atom.x, atom.y, atom.z});
                (*chainMap)[atom->chainID].push_back(newResidue);
            }
            delete atom;
        }
    }

    pdbFile.close();

    return *chainMap;
}



void saveToPDB(const fs::path& outputFilename, const std::map<char, std::vector<Residue*>>& chainMap) {
    std::ofstream pdbFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }

    // Set the formatting for residue.resSeq
    pdbFile << std::fixed << std::setprecision(0);

    // Write the modified atom records
    for (const auto& chainEntry : chainMap) {
        for (const Residue* residue : chainEntry.second) {
            for (const Atom& atom : residue->atoms) {
                pdbFile << "ATOM  ";
                pdbFile.width(5);
                pdbFile << std::right << atom.serial;
                pdbFile << "  ";
                pdbFile << std::setw(3) << std::left <<  atom.name;
                pdbFile << atom.altLoc << atom.resName;
                pdbFile << " ";
                pdbFile << chainEntry.first;
                pdbFile.width(4);
                pdbFile << std::right << residue->resSeq;
                //pdbFile << atom.iCode;
                pdbFile << "    ";
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[0];
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[1];
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[2];
                pdbFile.width(6);
                pdbFile.precision(2);
                pdbFile << std::fixed << atom.occupancy;
                pdbFile.width(6);
                pdbFile.precision(2);
                pdbFile << std::fixed << atom.tempFactor;
                pdbFile << "          ";
                pdbFile << atom.element << std::endl;
            }
        }
    }

    pdbFile.close();
}


void saveToPDBWithComments(const fs::path& outputFilename, const std::map<char, std::vector<Residue*>>& chainMap, std::vector<std::string>& comments) {
    std::ofstream pdbFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }

    // Set the formatting for residue.resSeq
    pdbFile << std::fixed << std::setprecision(0);

    for(const std::string& comment: comments)
    {
        pdbFile << "REMARK ";
        pdbFile.width(3);
        pdbFile << 999 << " ";
        pdbFile << comment << std::endl;
    }
    // Write the modified atom records
    for (const auto& chainEntry : chainMap) {
        for (const Residue* residue : chainEntry.second) {
            for (const Atom& atom : residue->atoms) {
                pdbFile << "ATOM  ";
                pdbFile.width(5);
                pdbFile << std::right << atom.serial;
                if(atom.name.size() == 4)
                {
                    pdbFile << " ";
                    pdbFile << std::setw(4) << std::left << atom.name;
                }
                else {
                    pdbFile << "  ";
                    pdbFile << std::setw(3) << std::left << atom.name;
                }
                pdbFile << atom.altLoc << atom.resName;
                pdbFile << " ";
                pdbFile << chainEntry.first;
                pdbFile.width(4);
                pdbFile << std::right << residue->resSeq;
                //pdbFile << atom.iCode;
                pdbFile << "    ";
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[0];
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[1];
                pdbFile.width(8);
                pdbFile.precision(3);
                pdbFile << std::fixed << atom.coords[2];
                pdbFile.width(6);
                pdbFile.precision(2);
                pdbFile << std::fixed << atom.occupancy;
                pdbFile.width(6);
                pdbFile.precision(2);
                pdbFile << std::fixed << atom.tempFactor;
                pdbFile << "          ";
                pdbFile << atom.element << std::endl;
            }
        }
    }

    pdbFile.close();
}
