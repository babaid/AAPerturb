//
// Created by babaid on 09.09.23.
//

#ifndef AA_SIDECHAIN_PERTURBATION_CPP_MONTECARLO_H
#define AA_SIDECHAIN_PERTURBATION_CPP_MONTECARLO_H

#include<vector>
#include<map>
#include<utility>
#include "molecules.h"

std::pair<char ,std::vector<std::size_t>> chooseRandomInterfaceResidue(std::map<char, std::vector<Residue*>>& chainMap, const std::map<char, std::vector<int>>&);
void rotateResidueSidechainRandomly(std::map<char, std::vector<Residue*>>& , char , std::size_t );

#endif //AA_SIDECHAIN_PERTURBATION_CPP_MONTECARLO_H
