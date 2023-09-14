//
// Created by babaid on 14.09.23.
//
#include<string>
#include<iostream>
#include<filesystem>
#include<vector>
#include<iterator>
#include<algorithm>
#include<thread>
#include<argparse/argparse.hpp>
#include "../include/pdbparser.h"
#include "../include/io.h"

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("pdbcleaner");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output PDB which are perturbed randomly.");


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
        if(!fs::is_directory(input_dir))
        {
            std::cout << "The directory provided does not exist. Exiting..." << std::endl;
            std::exit(1);
        }
    }
    if (auto ofn = program.present("-o"))
    {
        output_dir = *ofn;
        if(!fs::is_directory(output_dir))
        {
            std::cout << "The directory provided does not exist. Creating directory " << output_dir << std::endl;
            fs::create_directory(output_dir);
        }
    }


    auto files = findInputFiles(input_dir);
    for (auto const& file:files) {
        auto clean_structure = parsePDBToBeCleaned(file);
        std::vector<std::string> comments;
        comments.push_back("This file was previously reindexed and the waters and the hydrogens were removed");
        saveToPDBWithComments(output_dir/file.filename(), clean_structure,  comments);
    }
    return 0;
}



