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

std::unique_ptr<std::map<char, std::vector<Residue>>> parsePDB(const  fs::path& filename, bool excludewaters, bool deprotonate) {

    //work in progress
    std::unique_ptr<std::map<char, std::vector<Residue>>> chainMap = std::make_unique<std::map<char, std::vector<Residue>>>();

    std::ifstream pdbFile(filename);

    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return chainMap;
    }

    std::string line;

    //long long int prevResSeq{-1}, residueCounter{-1};


    while (std::getline(pdbFile, line)) {
        if (line.compare(0, 4, "ATOM") == 0 || ((line.compare(0, 6, "HETATM") == 0) && !excludewaters) )  {
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
            if (chainMap->find(atom.chainID) == chainMap->end()) {
                (*chainMap)[atom.chainID] = std::vector<Residue>();
                std::cout << "There where " << atom.resSeq << " residues in chain " << atom.chainID << std::endl;
                //prevResSeq = 0;
            }


            // Check if this residue is already in the chain's residues
            bool found = false;
            for (auto& residue : chainMap->at(atom.chainID)) {
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
                //newResidue.atom_coords.push_back({atom.x, atom.y, atom.z});
                chainMap->at(atom.chainID).emplace_back(std::move(newResidue));
                //(*chainMap)[atom.chainID].push_back(newResidue);
            }
        }
    }

    pdbFile.close();

    return chainMap;
}

std::unique_ptr<std::map<char, std::vector<Residue>>>  parsePDBToBeCleaned(const  fs::path& filename, bool excludewaters, bool deprotonate) {

    //work in progress
    std::unique_ptr<std::map<char, std::vector<Residue>>> chainMap = std::make_unique<std::map<char, std::vector<Residue>>>();

    std::ifstream pdbFile(filename);

    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return chainMap;
    }

    std::string line;

    long long int prevResSeq{-1}, residueCounter{-1};
    long long int atomcntr{0};
    bool parsingFirstModel = false;
    bool parsingAtoms = false;

    while (std::getline(pdbFile, line)) {
        if (line.compare(0, 5, "MODEL") == 0 && !parsingFirstModel) {
            parsingAtoms = true; // Start parsing atoms when a "MODEL" is encountered
            parsingFirstModel = true;
            continue;
        }
        if (line.compare(0, 6, "ENDMDL") == 0 && parsingAtoms) {
            parsingAtoms = false; // End parsing atoms when "ENDMDL" is encountered
            continue;
        }
        if (parsingAtoms && (line.compare(0, 4, "ATOM") == 0 || ((line.compare(0, 6, "HETATM") == 0) && !excludewaters)) ) {
            Atom atom;
            //atom.serial = std::stoi(line.substr(6, 5));
            atom.name = line.substr(12, 4);
            atom.altLoc = line[16];
            atom.resName = line.substr(17, 3);
            atom.chainID = line[21];
            atom.resSeq = std::stoi(line.substr(22, 4))-1;
            atom.coords = {std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                            std::stod(line.substr(46, 8))};
            atom.occupancy = std::stod(line.substr(54, 6));
            atom.tempFactor = std::stod(line.substr(60, 6));
            atom.element = line.substr(76, 2);

            // Remove whitespace from the atom name and element
            atom.name.erase(std::remove_if(atom.name.begin(), atom.name.end(), ::isspace), atom.name.end());
            atom.element.erase(std::remove_if(atom.element.begin(), atom.element.end(), ::isspace), atom.element.end());

            //If deprotonate then we skip
            if ((atom.element == "H" && deprotonate)) continue;

            //My head hurts thinking about how many problems this cause me in the previous months.

            atom.serial = ++atomcntr;

            // Check if this chain is already in the map
            if (chainMap->find(atom.chainID) == chainMap->end()) {
                std::cout << "There where " << residueCounter << " residues in chain " << atom.chainID << std::endl;
                (*chainMap)[atom.chainID] = std::vector<Residue>();
                residueCounter = -1;
            }
            if (atom.resSeq != prevResSeq) {
                residueCounter++;
            }
            prevResSeq = atom.resSeq;

            // Check if this residue is already in the chain's residues
            bool found = false;
            for (auto &residue: chainMap->at(atom.chainID)) {
                if (residue.resSeq == atom.resSeq && residue.resName == atom.resName) {
                    residue.atoms.push_back(std::move(atom));
                    found = true;
                    break;
                }
            }

            // If the residue doesn't exist, create a new one
            if (!found) {
                Residue newResidue;
                newResidue.chainID = atom.chainID;
                newResidue.resSeq = residueCounter;
                newResidue.resName = atom.resName;
                newResidue.atoms.emplace_back(std::move(atom));
                chainMap->at(atom.chainID).emplace_back(std::move(newResidue));

            }
        }

    }
    if (!parsingFirstModel) {
        pdbFile.clear(); // Reset the end-of-file flag
        pdbFile.seekg(0, std::ios::beg); // Rewind to the beginning of the file
        parsingAtoms = true;
        while (std::getline(pdbFile, line)) {
            if (parsingAtoms && (line.compare(0, 4, "ATOM") == 0 || ((line.compare(0, 6, "HETATM") == 0) && !excludewaters)) ) {
                Atom atom;
                //atom.serial = std::stoi(line.substr(6, 5));
                atom.name = line.substr(12, 4);
                atom.altLoc = line[16];
                atom.resName = line.substr(17, 3);
                atom.chainID = line[21];
                atom.resSeq = std::stoi(line.substr(22, 4)) - 1;
                atom.coords = {std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                               std::stod(line.substr(46, 8))};
                atom.occupancy = std::stod(line.substr(54, 6));
                atom.tempFactor = std::stod(line.substr(60, 6));
                atom.element = line.substr(76, 2);

                // Remove whitespace from the atom name and element
                atom.name.erase(std::remove_if(atom.name.begin(), atom.name.end(), ::isspace), atom.name.end());
                atom.element.erase(std::remove_if(atom.element.begin(), atom.element.end(), ::isspace), atom.element.end());

                //If deprotonate then we skip
                if ((atom.element == "H" && deprotonate)) continue;

                //My head hurts thinking about how many problems this cause me in the previous months.

                atom.serial = ++atomcntr;

                // Check if this chain is already in the map
                if (chainMap->find(atom.chainID) == chainMap->end()) {
                    (*chainMap)[atom.chainID] = std::vector<Residue>();
                    residueCounter = -1;
                }
                if (atom.resSeq != prevResSeq) {
                    residueCounter++;
                }
                prevResSeq = atom.resSeq;

                // Check if this residue is already in the chain's residues
                bool found = false;
                for (auto &residue: chainMap->at(atom.chainID)) {
                    if (residue.resSeq == atom.resSeq && residue.resName == atom.resName) {
                        residue.atoms.push_back(std::move(atom));
                        found = true;
                        break;
                    }
                }

                // If the residue doesn't exist, create a new one
                if (!found) {
                    Residue newResidue;
                    newResidue.chainID = atom.chainID;
                    newResidue.resSeq = residueCounter;
                    newResidue.resName = atom.resName;
                    newResidue.atoms.emplace_back(std::move(atom));
                    chainMap->at(atom.chainID).emplace_back(std::move(newResidue));

                }
            }

        }


    }
        pdbFile.close();


    return chainMap;
}



void saveToPDB(const fs::path& outputFilename, const std::unique_ptr<std::map<char, std::vector<Residue>>>& chainMap) {
    std::ofstream pdbFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }

    // Set the formatting for residue.resSeq
    pdbFile << std::fixed << std::setprecision(0);
    pdbFile << "MODEL        1" << std::endl;

    // Write the modified atom records
    for (const auto& chainEntry : *chainMap) {
        for (auto& residue : chainEntry.second) {
            for (const Atom& atom : residue.atoms) {
                pdbFile << "ATOM  ";
                pdbFile.width(5);
                pdbFile << std::right << atom.serial;
                pdbFile << "  ";
                pdbFile << std::setw(3) << std::left <<  atom.name;
                pdbFile << atom.altLoc << atom.resName;
                pdbFile << " ";
                pdbFile << chainEntry.first;
                pdbFile.width(4);
                pdbFile << std::right << residue.resSeq;
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
    pdbFile << "ENDMDL"<< std::endl;
    pdbFile.close();
}


void saveToPDBWithComments(const fs::path& outputFilename, const std::unique_ptr<std::map<char, std::vector<Residue>>>& chainMap, std::vector<std::string>& comments) {
    std::ofstream pdbFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }

    // Set the formatting for residue.resSeq
    pdbFile << std::fixed << std::setprecision(0);
    pdbFile << "MODEL        1" << std::endl;
    for(const std::string& comment: comments)
    {
        pdbFile << "REMARK ";
        pdbFile.width(3);
        pdbFile << 999 << " ";
        pdbFile << comment << std::endl;
    }
    // Write the modified atom records
    for (const auto& chainEntry : *chainMap) {
        for (auto& residue : chainEntry.second) {
            for (const Atom& atom : residue.atoms) {
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
                pdbFile << std::right << residue.resSeq+1;
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
    pdbFile << "ENDMDL"<< std::endl;
    pdbFile.close();
}
