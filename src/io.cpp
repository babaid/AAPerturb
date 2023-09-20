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

std::vector<fs::path> findInputFiles(const fs::path& path, const std::string extension)
{
    std::vector<fs::path> input_files;
    for (auto const& file : fs::directory_iterator{fs::absolute(path)})
    {
        if(fs::path(file).extension() == extension)
        {
            input_files.emplace_back(file);
        }
    }
    return input_files;
}

