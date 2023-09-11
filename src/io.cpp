//
// Created by babaid on 09.09.23.
//
#include<iostream>
#include<vector>
#include<string>
#include<filesystem>
#include<functional>
#include<iterator>
#include "../include/io.h"

namespace fs = std::filesystem;

std::vector<fs::path> createFileBatches(const fs::path& path, const std::size_t& batch_size)
{
    std::vector<fs::path> input_files;
    for (auto const& file : fs::directory_iterator{path})
    {
        if(fs::path(file).extension() == ".pdb")
        {
            input_files.emplace_back(file);
        }
    }
    return input_files;
}

