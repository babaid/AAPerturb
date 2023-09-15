//Header for parsing PDB files

#ifndef AAPERTURB_PDBPARSER_H
#define AAPERTURB_PDBPARSER_H
#include<vector>
#include<string>
#include<map>
#include<filesystem>
#include<valarray>
#include "molecules.h"

namespace fs = std::filesystem;

std::unique_ptr<std::map<char, std::vector<Residue>>>  parsePDB(const fs::path& filename, bool=true, bool=true);
std::unique_ptr<std::map<char, std::vector<Residue>>> parsePDBToBeCleaned(const  fs::path& filename, bool excludewaters=true, bool deprotonate=true);
void  saveToPDB(const fs::path& filename, const std::unique_ptr<std::map<char, std::vector<Residue>>>&);
void  saveToPDBWithComments(const fs::path& filename, const std::unique_ptr<std::map<char, std::vector<Residue>>>&, std::vector<std::string>& comments);


#endif
