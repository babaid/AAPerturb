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

/*
 * This function chooses a random Residue, given a map of possible chains with possible residue indices.
 */
std::pair<char ,std::vector<std::size_t>> chooseRandomResidue(const std::map<char, std::vector<int>>& interface_residue_indices)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<> dist(0, interface_residue_indices.size()-1);
    auto chain = interface_residue_indices.begin();
    std::advance(chain,  dist(rng));
    //Choose random residue on chain
    //std::size_t nelements = 1; //we want to change one residue for now
    std::random_device dev2;
    std::mt19937 rng2(dev2());
    std::vector<std::size_t> residues;
    //this somehow kills the program so I will use something else
    //std::sample(chain->second.begin(), chain->second.end(), std::back_inserter(residues), nelements, std::mt19937(std::random_device{}()));
    std::uniform_int_distribution<> resdist(0, chain->second.size()-2);
    residues.emplace_back(chain->second.at(resdist(rng2)));
    return std::make_pair(chain->first, residues);

}

double rotateResidueSidechainRandomly(std::unique_ptr<std::map<char, std::vector<Residue>>> & structure, char chain, std::size_t resNum, bool verbose)
{
    if (verbose)
    {
        std::cout <<  std::endl << "Size of chain: " <<structure->at(chain).size() << std::endl;
    }
    std::random_device dev;
    std::mt19937 rng(dev());
    double angles = 10; // keep it small or change the clash cutoff, if not changed there could still be clashes...
                        // rule of thumb <10 clash_cutoff -> 0.21 (approx. hydrogen covalent radius)
                        // the greater the angle range gets, the greater should be the clash cutoff
                        // optionally we could differentiate between types of atoms at clashes, but is it worth it?

    std::uniform_real_distribution<double> dist( -angles, angles);
    std::string resName = structure->at(chain).at(resNum).resName;
    Residue ref_res(structure->at(chain).at(resNum));
    //Residue ref_res(structure->at(chain).at(resNum));
    //std::cout << "Changing residue: "<< chain << "/" <<resName<<":"<< resNum+1<< std::endl;
    double rmsd{0};
    unsigned int patience = 0;
    //we do not prefer GLY PRO and ALA. Maybe only for displacement
    if(resName!="GLY" && resName!= "PRO" && resName!= "ALA") {
        while (rmsd == 0 && patience < 3) {
            for (const std::string &axis: amino_acids::axes::AMINO_MAP.at(resName)) {
                auto it_substructure = std::find(amino_acids::atoms::AMINO_MAP.at(resName).begin(),
                                                 amino_acids::atoms::AMINO_MAP.at(resName).end(), axis);
                std::size_t index = std::distance(amino_acids::atoms::AMINO_MAP.at(resName).begin(), it_substructure);
                auto first = amino_acids::atoms::AMINO_MAP.at(resName).begin() + index;

                std::vector<std::string> sub_atoms = std::vector<std::string>(first, amino_acids::atoms::AMINO_MAP.at(
                        resName).end());
                std::valarray<double> rot_coords = findRotationAxis(structure->at(chain).at(resNum), axis);
                double vec_norm = std::sqrt(std::pow(rot_coords, 2).sum());
                double angle = dist(rng);

                for (Atom& atom: structure->at(chain).at(resNum).atoms)
                    if (std::count(sub_atoms.begin(), sub_atoms.end(), atom.name)) {
                        if (verbose) std::cout << "Performing a rotation after " << *first << std::endl;
                        rotateCoordinatesAroundAxis(atom.coords, rot_coords / vec_norm, angle);

                    }
            }

                std::vector<std::vector<double>> distance_matrix =
                        calculateLocalDistanceMatrix(structure, structure->at(chain).at(resNum));

                if (detect_clashes(distance_matrix, 0.21)) {
                    if (verbose) std::cout << "Atoms clashed, retrying..." << std::endl;
                    structure->at(chain).at(resNum) = ref_res;
                    patience++;
                    angles = angles / 10;
                    std::uniform_real_distribution<double> dist(-angles, angles);
                    continue;
                } else patience = 0;
                rmsd = calculateRMSD(ref_res.atoms, structure->at(chain).at(resNum).atoms);
        }
    }
    return rmsd;
}
