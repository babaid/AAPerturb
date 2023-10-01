//
// Created by babaid on 01.10.23.
//
#include<argparse/argparse.hpp>
#include "pdbparser.h"
#include<filesystem>
#include "benchmarks.h"

int main(int argc, char *argv[]){
    argparse::ArgumentParser program("run-benchmarks");
    program.add_argument("-v", "--verbose")
            .help("Enable verbose mode. It is useful if you want to know if something is happening.")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files for the benchmarks.");

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
    std::cout << "Directory: " << input_dir << std::endl;

    benchmark_pdbparser(input_dir, 10);

    return 0;
}