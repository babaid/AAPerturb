//
// Created by babaid on 12.09.23.
//
// Finds and makes lists of the interfaces of a protein or a bunch of proteins.

#include<string>
#include<iostream>
#include<filesystem>
#include<array>
#include<vector>
#include<random>
#include<iterator>
#include<algorithm>
#include<thread>
#include<functional>
#include<format>
#include<fstream>
#include<argparse/argparse.hpp>
#include "../include/pdbparser.h"
#include "../include/geometry.h"
#include "../include/montecarlo.h"
#include "../include/molecules.h"
#include "../include/io.h"
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("iffind");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to determine the interfaces of.");
    program.add_argument("-o", "--output-dir")
            .help("The directory containing the output .sel files which contain the chain and the number of the interface residue.");

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    fs::path input_dir, output_dir;
    if (auto ifn = program.present("-i"))
    {
        input_dir = *ifn;
    }
    if (auto ofn = program.present("-o"))
    {
        output_dir = *ofn;
    }


    if(fs::is_directory(input_dir) && fs::is_directory(output_dir)) {
        std::vector<fs::path> files = findInputFiles(input_dir);
        std::vector<std::map<char, std::vector<int>>> all;
        for (auto file: files) {
            std::map<char, std::vector<Residue *>> structure = parsePDB(file);
            const std::map<char, std::vector<int>> interface_residues = findInterfaceResidues(structure, 12.0);

            fs::path outputFilename = output_dir / file.replace_extension(".sel").filename();

            std::ofstream selFile(outputFilename);

            if (!selFile.is_open()) {
                std::cerr << "Error: Unable to open file " << outputFilename << std::endl;
                return 1;
            }
            for (const auto &chain: interface_residues) {
                for (const auto &resnum: interface_residues.at(chain.first)) {
                    selFile << "/" << chain.first << ":" << resnum << " ";
                }
                selFile << std::endl;
            }
            selFile.close();
        }
    }
    else if (fs::exists(input_dir) && !fs::is_directory(input_dir) && output_dir.empty())
    {
        std::map<char, std::vector<Residue *>> structure = parsePDB(input_dir);
        const std::map<char, std::vector<int>> interface_residues = findInterfaceResidues(structure, 12.0);
        for (const auto &chain: interface_residues) {
            for (const auto &resnum: interface_residues.at(chain.first)) {
                std::cout << "/" << chain.first << ":" << resnum << " ";
            }
            std::cout << std::endl;
        }
    }
    else
    {
        std::cerr << "Some error occurred. Try again." << std::endl;
        std::exit(1);
    }
    return 0;
}
