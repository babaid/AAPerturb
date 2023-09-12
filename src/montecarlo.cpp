//
// Created by babaid on 09.09.23.
//
#include<vector>
#include<map>
#include<random>
#include<algorithm>
#include<iostream>
#include<utility>

#include "../include/geometry.h"
#include "../include/montecarlo.h"
#include "../include/constants.h"


std::pair<char ,std::vector<std::size_t>> chooseRandomInterfaceResidue(std::map<char, std::vector<Residue*>>& chainMap, const std::map<char, std::vector<int>>& interface_residue_indices)
{
    //Find interface residues
    //std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(chainMap);
    //Choose a random chain
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::size_t> dist(0, chainMap.size()-1);
    auto chain = interface_residue_indices.begin();
    std::advance(chain,  dist(rng));
    //Choose random residue on chain
    //std::size_t nelements = 1; //we want to change one residue for now
    std::vector<std::size_t> residues;
    //this somehow kills the program so I will
    //std::sample(chain->second.begin(), chain->second.end(), std::back_inserter(residues), nelements, std::mt19937(std::random_device{}()));
    std::uniform_int_distribution<std::size_t> resdist(0, interface_residue_indices.size()-1);
    residues.emplace_back(resdist(rng));
    return std::make_pair(chain->first, residues);

}

void rotateResidueSidechainRandomly(std::map<char, std::vector<Residue*>>& structure, char chain, std::size_t resNum)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    double angles = 90;
    std::uniform_real_distribution<double> dist( -angles, angles);
    std::string resName = structure.at(chain).at(resNum)->atoms[0].resName;
    Residue* ref_res = new Residue();
    *ref_res = *structure.at(chain).at(resNum);
    //Residue ref_res(structure->at(chain).at(resNum));
    std::cout << "Changing residue: "<< chain << "/" <<resName<<":"<< resNum+1<< std::endl;
    double rmsd{0};
    unsigned int patience = 0;
    //we do not prefer GLY PRO and ALA. Maybe only for displacement
    if(resName!="GLY" && resName!= "PRO" && resName!= "ALA") {
        while (rmsd == 0 || patience > 3) {
            for (const std::string &axis: amino_acids::axes::AMINO_MAP.at(resName)) {
                std::cout << "Rotating around axis: " << axis << std::endl;
                auto it_substructure = std::find(amino_acids::atoms::AMINO_MAP.at(resName).begin(),
                                                 amino_acids::atoms::AMINO_MAP.at(resName).end(), axis);
                std::size_t index = std::distance(amino_acids::atoms::AMINO_MAP.at(resName).begin(), it_substructure);
                auto first = amino_acids::atoms::AMINO_MAP.at(resName).begin() + index;
                std::vector<std::string> sub_atoms = std::vector<std::string>(first, amino_acids::atoms::AMINO_MAP.at(
                        resName).end());

                std::cout << "Rotating atoms" << std::endl;
                for (const std::string &s: sub_atoms) std::cout << s << " ";
                std::cout << std::endl;

                std::valarray<double> rot_coords = findRotationAxis(structure.at(chain).at(resNum), axis);

                double vec_norm = std::sqrt(std::pow(rot_coords, 2).sum());
                double angle = dist(rng);

                for (Atom& atom: structure[chain][resNum]->atoms)
                    if (std::count(sub_atoms.begin(), sub_atoms.end(), atom.name)) {

                        rotateCoordinatesAroundAxis(atom.coords, rot_coords / vec_norm, angle);

                    }
            }
            auto distance_matrix = calculateLocalDistanceMatrix(structure, structure.at(chain).at(resNum));
           
            if (detect_clashes(distance_matrix, 0.21))  {
                for (auto row:distance_matrix)
                {
                    for(auto el:row) std::cout<< el << " ";
                    std::cout << std::endl;
                }
                std::cout << "Atoms clashed, retrying..." << std::endl;
                *structure.at(chain).at(resNum) = *ref_res;
                patience++;
                break;
            }
            patience=0;
            rmsd = calculateRMSD(ref_res->atoms, structure.at(chain).at(resNum)->atoms);
            std::cout << "Current RMSD of the residue is: " << rmsd << std::endl;
        }
    }
    delete ref_res;
}
