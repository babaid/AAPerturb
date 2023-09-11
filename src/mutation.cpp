//
// Created by babaid on 09.09.23.
//


#include "../include/mutation.h"
#include "../include/constants.h"
#include "../include/geometry.h"
#include<algorithm>
#include<iostream>
#include<cmath>
void mutateResidue(std::map<char, std::vector<Residue>>& chains, char chain, std::size_t resId, std::string to){
    --resId;
    std::string from = chains.at(chain).at(resId).resName;
    if(from == to)
    {
        std::cout << "This muattion will have no effect, as the two residues are the same." << std:: endl;
    }
    if (to == "GLY")
    {
        auto endIt = std::remove_if(
                chains.at(chain).at(resId).atoms.begin(),
                chains.at(chain).at(resId).atoms.end(),
                [&to=to](const Atom& atom){
                    return std::find(amino_acids::atoms::AMINO_MAP.at(to).begin(),
                                     amino_acids::atoms::AMINO_MAP.at(to).end(),
                                     atom.name) == amino_acids::atoms::AMINO_MAP.at(to).end();
                });
        chains.at(chain).at(resId).atoms.erase(endIt, chains.at(chain).at(resId).atoms.end());
        //renaming residue
        chains.at(chain).at(resId).resName = "GLY";
        for (Atom& atom:chains.at(chain).at(resId).atoms)
        {
            atom.resName = "GLY";
        }

    }
    else if (to == "ALA")
    {
        if (from != "GLY") {
            auto endIt = std::remove_if(
                    chains.at(chain).at(resId).atoms.begin(),
                    chains.at(chain).at(resId).atoms.end(),
                    [&to = to](const Atom &atom) {
                        return std::find(amino_acids::atoms::AMINO_MAP.at(to).begin(),
                                         amino_acids::atoms::AMINO_MAP.at(to).end(),
                                         atom.name) == amino_acids::atoms::AMINO_MAP.at(to).end();
                    });
            chains.at(chain).at(resId).atoms.erase(endIt, chains.at(chain).at(resId).atoms.end());
            //renaming residue
            chains.at(chain).at(resId).resName = "ALA";
            for (Atom &atom: chains.at(chain).at(resId).atoms) {
                atom.resName = "ALA";
            }
        }
        else
        {
            //add a CB atom at 1.525 AA over CA atom in the N-CA-C plane, rotate it into having a 111.1 angle with N-CA-CB and a 109 angle with C-CA-CB
            std::vector<Atom>::const_iterator n_it = std::find_if(chains.at(chain).at(resId).atoms.begin(),
                                                                chains.at(chain).at(resId).atoms.end(),
                                                                [](const Atom& atom)->bool {return (atom.name == "N");} );
            std::vector<Atom>::const_iterator ca_it = std::find_if(chains.at(chain).at(resId).atoms.begin(),
                                                                   chains.at(chain).at(resId).atoms.end(),
                                                                   [](const Atom& atom)->bool {return (atom.name == "CA");} );
            std::vector<Atom>::const_iterator c_it = std::find_if(chains.at(chain).at(resId).atoms.begin(),
                                                                   chains.at(chain).at(resId).atoms.end(),
                                                                   [](const Atom& atom)->bool {return (atom.name == "C");} );
            std::valarray<double> perp_vec = crossProduct(n_it->coords - ca_it->coords, c_it->coords-ca_it->coords);
            perp_vec = perp_vec/norm(perp_vec);
            Atom CB{*ca_it};
            CB.name = "CB";
            //move to the right distance
            CB.coords += c_it->coords+n_it->coords;

            //rotate into place

            //std::valarray<double> n_c = (n_it->coords+2*ca_it->coords);
            //n_c /= norm(n_c);
            //std::valarray<double> c = (c_it->coords+2*ca_it->coords);
            //c /= norm(c);

            //rotateCoordinatesAroundAxis(CB.coords, n_c, -90);
            //rotateCoordinatesAroundAxis(CB.coords, c, -90);

            //CB.coords += perp_vec;
            //insert the thing after the CA atom
            chains.at(chain).at(resId).atoms.insert(ca_it, CB);
            //need to renumber every atom.
            renumberAtoms(chains);

        }
    }
    //continue wth more
}


void renumberAtoms(std::map<char, std::vector<Residue>>& chainMap) {
    int new_serial = 1;
    for (auto &chainEntry: chainMap) {
        for (Residue &residue: chainEntry.second) {
            for (Atom &atom: residue.atoms) {
                atom.serial = new_serial;
                ++new_serial;
            }
        }
    }
}
