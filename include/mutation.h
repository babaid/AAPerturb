//Prototype of mutations
// Created by babaid on 09.09.23.
//

#ifndef AA_SIDECHAIN_PERTURBATION_CPP_MUTATION_H
#define AA_SIDECHAIN_PERTURBATION_CPP_MUTATION_H

#include<map>
#include<vector>
#include<string>
#include "pdbparser.h"

void mutateResidue(std::map<char, std::vector<Residue>>& chains, char chain, std::size_t resId, std::string to);
void renumberAtoms(std::map<char, std::vector<Residue>>& chainMap);


#endif //AA_SIDECHAIN_PERTURBATION_CPP_MUTATION_H
