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
std::pair<char , std::size_t> chooseRandomResidue(const std::map<char, std::vector<int>>& interface_residue_indices)
{
    //This thread local stuff is something really cool
    // Need to explain it here, otherwise I forget it.
    // If using RNGs with the same seed, random out of range errors can occur, as they share their state amongst threads
    // Using the thread_local keyword they become independent, do not share the state anymore and there is no synchronization required.
    // An other option is to use global RNGs, but this is a bad idea because race conditions will probably occur,
    // and you have to use mutexes to synchronize the r.n. generation.
    thread_local std::random_device thread_dev;
    thread_local std::mt19937 thread_rng(thread_dev());
    std::size_t resindex = std::numeric_limits<std::size_t>::max();

    std::uniform_int_distribution<> dist(0, interface_residue_indices.size() - 1); // for chain selection
    auto chain = interface_residue_indices.begin();
    if (!interface_residue_indices.empty()) {
        while(resindex == std::numeric_limits<std::size_t>::max()) {
            chain = interface_residue_indices.begin();
            std::advance(chain, dist(thread_rng));
            if (chain != interface_residue_indices.end() && !chain->second.empty()) {
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

double rotateResidueSidechainRandomly(std::unique_ptr<std::map<char, std::vector<Residue>>> & structure,  char& chain, std::size_t resNum, bool verbose)
{

    if (verbose)
    {
        std::cout <<  std::endl << "Size of chain: " <<structure->at(chain).size() << std::endl << "Trying to perturb..." << std::endl;
    }

    thread_local std::random_device thread_dev;
    thread_local std::mt19937 thread_rng(thread_dev());

    double angles = 10; // keep it small or change the clash cutoff, if not changed there could still be clashes...
                        // rule of thumb <10 clash_cutoff -> 0.21 (approx. hydrogen covalent radius)
                        // the greater the angle range gets, the greater should be the clash cutoff
                        // optionally we could differentiate between types of atoms at clashes, but is it worth it?
    std::uniform_real_distribution<double> dist( -angles, angles);
    Residue ref_res(structure->at(chain).at(resNum));
    std::string resName = structure->at(chain).at(resNum).resName;
    ref_res = structure->at(chain).at(resNum);


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
                double angle = dist(thread_rng);

                for (Atom& atom: structure->at(chain).at(resNum).atoms)
                    if (std::count(sub_atoms.begin(), sub_atoms.end(), atom.name)) {
                        //if (verbose) std::cout << "Performing a rotation after " << *first << std::endl;
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
