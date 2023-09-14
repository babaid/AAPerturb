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

std::map<char, std::vector<std::unique_ptr<Residue>>>  parsePDB(const fs::path& filename, bool=true, bool=true);
std::map<char, std::vector<std::unique_ptr<Residue>>> parsePDBToBeCleaned(const  fs::path& filename, bool excludewaters=true, bool deprotonate=true);
void  saveToPDB(const fs::path& filename, const std::map<char, std::vector<std::unique_ptr<Residue>>>& chainMap);
void  saveToPDBWithComments(const fs::path& filename, const std::map<char, std::vector<std::unique_ptr<Residue>>>& chainMap, std::vector<std::string>& comments);


#endif
