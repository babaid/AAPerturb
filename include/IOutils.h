//Header for handling input and output files/creating directory structures
// Created by babaid on 09.09.23.
//

#ifndef AAPERTURB_IOUTILS_H
#define AAPERTURB_IOUTILS_H
#include<vector>
#include<filesystem>
namespace fs = std::filesystem;

std::vector<fs::path> findInputFiles(const fs::path& path, const std::string extension=  ".pdb");
void saveMatrixAsTSV(std::vector<std::vector<double>>&, fs::path);
std::size_t number_of_files_in_directory(fs::path path);
#endif //AA_SIDECHAIN_PERTURBATION_CPP_IO_H
