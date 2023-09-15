//
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB__MONTECARLO_H
#define AAPERTURB__MONTECARLO_H

#include<vector>
#include<map>
#include<utility>
#include "molecules.h"

std::pair<char , std::size_t> chooseRandomResidue(const std::map<char, std::vector<int>>& interface_residue_indices);
double rotateResidueSidechainRandomly(std::unique_ptr<std::map<char, std::vector<Residue>>> & , char , std::size_t, bool);

#endif
