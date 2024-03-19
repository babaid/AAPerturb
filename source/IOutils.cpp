#include "IOutils.h"

#include<iostream>
#include<vector>
#include<string>
#include<filesystem>
#include<functional>
#include<iterator>
#include<fstream>


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

void saveMatrixAsTSV(std::vector<std::vector<double>> & mat, fs::path outputFilename) {
    std::ofstream TSVFile(outputFilename);
    if (!TSVFile.is_open()) {
        std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
        return;
    }
    // Set the formatting for residue.resSeq
    TSVFile << std::fixed << std::setprecision(2);
    for(auto const& row: mat)
    {
        for(auto const& element:row)
        {
            TSVFile << element << "\t";
        }
        TSVFile << std::endl;
    }
    TSVFile.close();
}

std::size_t number_of_files_in_directory(fs::path path)
{
    using std::filesystem::directory_iterator;
    using fp = bool (*)( const std::filesystem::path&);
    return std::count_if(directory_iterator(path), directory_iterator{}, (fp)std::filesystem::is_regular_file);
}

