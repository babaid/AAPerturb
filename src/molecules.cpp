//
// Created by babaid on 09.09.23.
//
#include<iostream>
#include<algorithm>
#include<string>
#include<fstream>
#include<filesystem>
#include<exception>
#include<random>
#include "io.h"
#include "molecules.h"
#include "geometry.h"
#include "constants.h"

namespace fs=std::filesystem;


bool Atom::operator==(const Atom &atom) const {
    bool ser = this->serial == atom.serial;
    bool name = this->name == atom.name;
    bool rname = this->resName == atom.resName;
    bool el = this->element == atom.element;
    return ser && name && rname && el;
}

double calculateDistance(const Atom& atom1, const Atom& atom2) {
    return std::sqrt(std::pow(atom2.coords - atom1.coords, 2).sum());
}


PDBStructure::PDBStructure(const fs::path filename, bool verbose): verbose(verbose) {
    //work in progress
    std::ifstream pdbFile(filename);

    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        throw;
    }
    std::string line;
    while (std::getline(pdbFile, line)) {
        if (line.compare(0, 4, "ATOM") == 0)  {
            numAtoms++;
            Atom atom;
            atom.serial = std::stoi(line.substr(6, 5));
            atom.name = line.substr(12, 4);
            atom.altLoc = line[16];
            atom.resName = line.substr(17, 3);
            atom.chainID = line[21];
            if(!atom.chainID) throw;
            atom.resSeq = std::stoi(line.substr(22, 4))-1;
            try {
                atom.coords = {std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                               std::stod(line.substr(46, 8))};
                atom.occupancy = std::stod(line.substr(54, 6));
                atom.tempFactor = std::stod(line.substr(60, 6));
                atom.element = line.substr(76, 4);
            }
            catch (...)
            {

                std::cerr << "Something is obscured in the PDB file. Maybe run the cleaner or check if there is anything wrong." << std::endl;
                throw;
            }
            // Remove whitespace from the atom name
            atom.element.erase(std::remove_if(atom.element.begin(), atom.element.end(), ::isspace), atom.element.end());
            atom.name.erase(std::remove_if(atom.name.begin(), atom.name.end(), ::isspace), atom.name.end());

            // Check if this chain is already in the map
            if (chains.find(atom.chainID) == chains.end()) {
                chains[atom.chainID] = std::vector<Residue>();
            }

            // Check if this residue is already in the chain's residues
            bool found = false;
            for (auto& residue : chains.at(atom.chainID)) {
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
                chains.at(atom.chainID).emplace_back(std::move(newResidue));
                //(*chainMap)[atom.chainID].push_back(newResidue);
            }
        }
    }
    pdbFile.close();
}

void PDBStructure::calculateDistanceMatrix() {
    dist_mat =  std::vector<std::vector<double>>(numAtoms, std::vector<double>(numAtoms, 0.0));
    int currentIndex = 0;
    // Iterate through the chainMap structure
    for (const auto& chainEntry1 : chains) {
        for (auto const& residue1 : chainEntry1.second) {
            for (const Atom& atom1 : residue1.atoms) {
                int otherIndex = 0;
                for (const auto& chainEntry2 : chains) {
                    for (auto const&  residue2 : chainEntry2.second) {
                        for (const Atom& atom2 : residue2.atoms) {
                            double distance = calculateDistance(atom1, atom2);
                            dist_mat[currentIndex][otherIndex] = distance;
                            if (currentIndex == otherIndex) {
                                dist_mat.at(currentIndex).at(otherIndex) = std::numeric_limits<double>::max();
                            }
                            otherIndex++;
                        }
                    }
                }
                currentIndex++;
            }
        }
    }
}

std::vector<std::vector<double>> PDBStructure::calculateLocalDistanceMatrix(Residue &refres) {

    std::vector<std::vector<double>> distanceMatrix =  std::vector<std::vector<double>>(numAtoms, std::vector<double>(refres.atoms.size(), 0.0));

    int currentIndex = 0;

    // Iterate through the chainMap structure
    for (const auto& chainEntry1 : chains) {
        for (auto const& residue1: chainEntry1.second) {
            for (const Atom& atom1: residue1.atoms) {
                int otherIndex = 0;
                for (const Atom& atom2: refres.atoms) {
                    double distance = calculateDistance(atom1, atom2);
                    distanceMatrix.at(currentIndex).at(otherIndex) = distance;
                    // Set diagonal elements to max
                    if (residue1.resSeq == refres.resSeq && atom1 == atom2) {
                        distanceMatrix.at(currentIndex).at(otherIndex) = std::numeric_limits<double>::max();
                    }
                    otherIndex++;
                }
                currentIndex++;
            }
        }
    }
    return distanceMatrix;
}


void PDBStructure::updateDistanceMatrixLocally(std::vector<std::vector<double>>& newPart, char chain, unsigned int resNum){
    auto atompos = chains.at(chain).at(resNum).atoms.at(0).serial-1;
    for(std::size_t i{0}; i<numAtoms; i++)
    {
        for(int j{atompos}; j<atompos+newPart.at(0).size(); ++j){
            dist_mat.at(i).at(j) = newPart.at(i).at(j-atompos);
            dist_mat.at(j).at(i) = newPart.at(i).at(j-atompos); // Symmetry. Maybe not optimal. IDK
        }
    }
}

std::vector<double> one_hot(unsigned int index, unsigned int num_classes) {
    std::vector<double> vec(num_classes, 0.0);
    vec.at(index) = 1.0;
    return vec;
}


void PDBStructure::savePDBStructure(fs::path outputFilename, const std::vector<std::string>& comments) {
    std::ofstream pdbFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }pdbFile << "MODEL        1" << std::endl;
    // Set the formatting for residue.resSeq
    pdbFile << std::fixed << std::setprecision(0);
    //pdbFile << "MODEL        1" << std::endl;

    for(const std::string& comment: comments)
    {
        pdbFile << "REMARK ";
        pdbFile.width(3);
        pdbFile << 999 << " ";
        pdbFile << comment << std::endl;
    }
                pdbFile.width(5);
    // Write the modified atom records
    for (const auto& chainEntry : chains) {
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
    //pdbFile << "ENDMDL"<< std::endl;
    pdbFile.close();
}

void PDBStructure::calculateAtomicFeatureMatrix() {
    std::vector<std::vector<double>> feat_mat(numAtoms, std::vector<double>(12));
    for (const auto& chainEntry : chains) {
        for (auto &residue: chainEntry.second) {
            for (const Atom &atom: residue.atoms) {
                auto element = one_hot(atoms::ELEMENTS.at(atom.element), 8);
                element.push_back(atoms::ELECTRONEGATIVITIES.at(atom.element));
                element.push_back(atoms::N_VALENCE_ELECTRONS.at(atom.element));
                element.push_back(atoms::FIRST_IONIZATION_ENERGIES.at(atom.element));
                element.push_back(atoms::POLARIZABILITIES.at(atom.element));
                feat_mat.push_back(element);
            }
        }


    }
    atomic_feat_mat = std::move(feat_mat);

}

void PDBStructure::saveFeatureMat(fs::path outputPath, bool flush) {
    auto outputFilename = outputPath/"at_feat.tsv";
    if(!atomic_feat_mat.empty()) {
        saveMatrixAsTSV(atomic_feat_mat, outputFilename);
        if (flush) {
            atomic_feat_mat.clear();
        }
    } else {
        std::cerr << "Before saving a feature matrix, you should actually have one to begin with." << std::endl
        << "I will forget about this and move on like nothing happened, but beware." << std::endl;
    }
}

void PDBStructure::saveCoords(fs::path outputPath) {
    auto outputFilename = outputPath/"coords.tsv";
    std::ofstream coordinateFile(outputFilename);
    //std::ofstream pdbFile(outputFilename, std::ios::out | std::ios::binary);
    if (!coordinateFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }
    for (const auto& chainEntry : chains) {
        for (auto &residue: chainEntry.second) {
            for (const Atom &atom: residue.atoms) {
                for(auto const& coord: atom.coords)
                {
                    coordinateFile << coord << "\t";
                }
                coordinateFile << std::endl;
            }
        }


    }

}

void PDBStructure::getNumberOfResidues() {
    std::cout << "Residue counts in each chain: " << std::endl;
    auto print_chain_n = [](auto const& elem){std::cout << elem.first << ": " << elem.second.size() << ", ";};
    std::for_each(chains.begin(), chains.end(), print_chain_n);
    std::cout << std::endl;
}

void PDBStructure::findInterfaceResidues(double cutoff) {
    std::map<char, std::vector<unsigned>> ifres;
    for (const auto& chainEntry1 : chains) {
        ifres[chainEntry1.first] = std::vector<unsigned>();
        for (auto const& residue1 : chainEntry1.second) {
            for (const auto& chainEntry2 : chains) {
                if (chainEntry1.first != chainEntry2.first) {
                    for (auto const & residue2 : chainEntry2.second) {
                        if (areResiduesNeighbors(residue1, residue2, cutoff)) {
                            // Check if residues are adjacent by comparing residue sequence numbers
                            ifres[chainEntry1.first].emplace_back(residue1.resSeq-1);
                            //interfaceResidues.push_back(residue1.resSeq);
                            break;  // No need to check further for this residue1
                        }
                    }
                }
            }
        }
    }
    interfaceResidues = std::move(ifres);
}

void PDBStructure::getInterfaceResidues() {
    std::cout << "Found following number of interface residues in the chains: ";
    auto print_chain_n = [](auto const& elem){std::cout << elem.first << ": " << elem.second.size() << ", ";};
    std::for_each(interfaceResidues.begin(), interfaceResidues.end(), print_chain_n);
    std::cout<< std::endl;
}

Residue PDBStructure::getResidue(char chain, unsigned int resSeq) {
    return chains.at(chain).at(resSeq);
}

void PDBStructure::setResidue(const Residue & res) {
    chains.at(res.chainID).at(res.resSeq) = std::move(res);
}

std::vector<std::vector<double>> PDBStructure::getDistMat() {
    return dist_mat;
}

void PDBStructure::saveDistMat(fs::path outputFilename) {
    saveMatrixAsTSV(dist_mat, outputFilename);

}

std::pair<char, std::size_t> PDBStructure::chooseRandomResidue() const {
    thread_local std::random_device thread_dev;
    thread_local std::mt19937 thread_rng(thread_dev());
    std::size_t resindex = std::numeric_limits<std::size_t>::max();

    std::uniform_int_distribution<> dist(0, interfaceResidues.size() - 1); // for chain selection
    auto chain = interfaceResidues.begin();
    if (!interfaceResidues.empty()) {
        while(resindex == std::numeric_limits<std::size_t>::max()) {
            chain = interfaceResidues.begin();
            std::advance(chain, dist(thread_rng));
            if (chain != interfaceResidues.end() && !chain->second.empty()) {
                std::uniform_int_distribution<std::size_t> elementDist(0, chain->second.size() - 1);
                resindex = chain->second.at(elementDist(thread_rng));
            }
            else
            {
                throw std::out_of_range("Some issue at random residue selection, either chain is empty or there are no chains.");
            }

        }
        return std::make_pair(chain->first, resindex);
    }
    else
    {
        throw std::out_of_range("Something is terribly wrong with the residue interfaces you try to feed me. They seem to be empty. Maybe increase the cutoff.");
    }

}

double PDBStructure::rotateResidueSidechainRandomly(char chain, std::size_t resNum) {
    if (verbose)
    {
        std::cout <<  std::endl << "Size of chain: " << chains.at(chain).size() << std::endl << "Trying to perturb..." << std::endl;
    }

    thread_local std::random_device thread_dev;
    thread_local std::mt19937 thread_rng(thread_dev());

    //this is where you could use your own distribution of angles
    double angles = 10; // keep it small or change the clash cutoff, if not changed there could still be clashes...
    // rule of thumb <10 clash_cutoff -> 0.21 (approx. hydrogen covalent radius)
    // the greater the angle range gets, the greater should be the clash cutoff
    // optionally we could differentiate between types of atoms at clashes, but is it worth it?
    std::uniform_real_distribution<double> dist( -angles, angles);

    Residue ref_res(chains.at(chain).at(resNum));
    std::string resName = chains.at(chain).at(resNum).resName;
    ref_res = chains.at(chain).at(resNum);


    //std::cout << "Changing residue: "<< chain << "/" <<resName<<":"<< resNum+1<< std::endl;
    double rmsd{0};
    //we do not prefer GLY PRO and ALA. Maybe only for displacement
    if(resName!="GLY" && resName!= "PRO" && resName!= "ALA" && amino_acids::axes::AMINO_MAP.find(resName) != amino_acids::axes::AMINO_MAP.end() ) {
        for (const std::string &axis: amino_acids::axes::AMINO_MAP.at(resName)) {
            auto it_substructure = std::find(amino_acids::atoms::AMINO_MAP.at(resName).begin(),
                                             amino_acids::atoms::AMINO_MAP.at(resName).end(), axis);
            std::size_t index = std::distance(amino_acids::atoms::AMINO_MAP.at(resName).begin(), it_substructure);
            auto first = amino_acids::atoms::AMINO_MAP.at(resName).begin() + index;

            std::vector<std::string> sub_atoms = std::vector<std::string>(first, amino_acids::atoms::AMINO_MAP.at(
                    resName).end());
            std::valarray<double> rot_coords = findRotationAxis(chains.at(chain).at(resNum), axis);
            double vec_norm = std::sqrt(std::pow(rot_coords, 2).sum());
            double angle = dist(thread_rng);

            for (Atom& atom: chains.at(chain).at(resNum).atoms)
                if (std::count(sub_atoms.begin(), sub_atoms.end(), atom.name)) {
                    //if (verbose) std::cout << "Performing a rotation after " << *first << std::endl;
                    rotateCoordinatesAroundAxis(atom.coords, rot_coords / vec_norm, angle);

                }
        }

        auto local_distance_matrix = this->calculateLocalDistanceMatrix(chains.at(chain).at(resNum));

        if (detect_clashes(local_distance_matrix, 0.21)) {
            if (verbose) std::cout << "Atoms clashed, retrying..." << std::endl;
            chains.at(chain).at(resNum) = ref_res;
            angles = angles / 10;
            std::uniform_real_distribution<double> dist(-angles, angles);
        } else this->updateDistanceMatrixLocally(local_distance_matrix, chain, resNum); // update distance matrix

        rmsd = calculateRMSD(ref_res.atoms, chains.at(chain).at(resNum).atoms);

    }
    return rmsd;
}

void PDBStructure::saveLocalDistMat(fs::path outputFile, char chain, int resNum) {
    auto local_dist_mat = this->calculateLocalDistanceMatrix(chains.at(chain).at(resNum));
    saveMatrixAsTSV(local_dist_mat, outputFile);
}

