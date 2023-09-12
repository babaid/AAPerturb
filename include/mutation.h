//Prototype of mutations
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB_MUTATION_H
#define AAPERTURB__MUTATION_H

#include<map>
#include<vector>
#include<string>
#include "pdbparser.h"

void mutateResidue(std::map<char, std::vector<Residue>>& chains, char chain, std::size_t resId, std::string to);
void renumberAtoms(std::map<char, std::vector<Residue>>& chainMap);


#endif
