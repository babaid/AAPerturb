//Header for handling input and output files/creating directory structures
// Created by babaid on 09.09.23.
//

#ifndef AA_SIDECHAIN_PERTURBATION_CPP_IO_H
#define AA_SIDECHAIN_PERTURBATION_CPP_IO_H
#include<vector>
#include<filesystem>
namespace fs = std::filesystem;

std::vector<fs::path> createFileBatches(const fs::path& path, const std::size_t& batch_size);

#endif //AA_SIDECHAIN_PERTURBATION_CPP_IO_H
