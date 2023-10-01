//
// Created by babaid on 01.10.23.
//

#include<argparse/argparse.hpp>
#include "pdbparser.h"
#include "io.h"
#include<filesystem>
#include<chrono>
#include "benchmarks.h"
namespace fs = std::filesystem;

void benchmark_pdbparser(fs::path& input_path, unsigned cnt)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto files = findInputFiles(input_path);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time spent on finding all the input files: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start) << std::endl;
    std::chrono::milliseconds duration{0};
    for (unsigned i{0}; i<cnt; i++){
        for (auto& file:files) {
            start = std::chrono::high_resolution_clock::now();
            auto structure = parsePDB(file);
            end = std::chrono::high_resolution_clock::now();
            duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        }
        duration/=files.size();
    }
    std::cout << "Overall avg. time for parsing a PDB file: " << duration << std::endl;
}
