//Header for parsing PDB files

#ifndef PDBPARSER_H
#define PDBPARSER_H
#include<vector>
#include<string>
#include<map>
#include<valarray>
#include "molecules.h"




std::map<char, std::vector<Residue*>>  parsePDB(const std::string& filename, bool=true, bool=true);
void  saveToPDB(const std::string& filename, const std::map<char, std::vector<Residue*>>& chainMap);
bool checkParsedPDB(const std::map<char, std::vector<Residue*>>&);

#endif
