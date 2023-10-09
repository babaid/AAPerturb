//
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB__MONTECARLO_H
#define AAPERTURB__MONTECARLO_H

#include<vector>
#include<map>
#include<utility>
#include "molecules.h"

class PDBStructure;
struct Residue;

std::pair<char , std::size_t> chooseRandomResidue(const std::unique_ptr<PDBStructure>& structure);
double rotateResidueSidechainRandomly(std::unique_ptr<std::map<char, std::vector<Residue>>> & , char& , std::size_t, bool);

#endif
