//
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB__MONTECARLO_H
#define AAPERTURB__MONTECARLO_H

#include<vector>
#include<map>
#include<utility>
#include "molecules.h"

std::pair<char ,std::vector<std::size_t>> chooseRandomResidue(const std::map<char, std::vector<int>>&);
double rotateResidueSidechainRandomly(std::map<char, std::vector<std::unique_ptr<Residue>>>& , char , std::size_t );

#endif
