//Main executable...
// Created by babaid on 07.09.23.
//
#include "perturbrun.h"
#include "spdlog/spdlog.h"

#include<argparse/argparse.hpp>

#include<iostream>
#include<filesystem>
#include<iterator>
#include<cmath>



//Verbose mode is currently not thread safe. I need to use mutexes or something...
using namespace std::chrono_literals;
namespace fs = std::filesystem;


int main(int argc, char *argv[]) {

    argparse::ArgumentParser program("aaperturb", "2.0.0");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output PDB which are perturbed randomly.");
    program.add_argument("-b", "--batch-size")
            .scan<'d', std::size_t>()
            .default_value(std::size_t(1))
            .help("The number of batches for multiprocessing");
    program.add_argument("-N", "--num-variations")
            .scan<'d', std::size_t>()
            .default_value(std::size_t(1))
            .help("The number of variations of a single protein.");

    program.add_argument("-q", "--max-bbangle")
            .scan<'g', double>()
            .default_value(double(0.5))
            .help("The maximal backbone rotation angle.");

    program.add_argument("-w", "--max-schangle")
            .scan<'g', double>()
            .default_value(double(0.5))
            .help("The maximal sidechain rotation angle.");


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
            spdlog::error("The directory provided does not exist. Exiting...");
            std::exit(1);
        }
    }
    if (auto ofn = program.present("-o"))
    {
        output_dir = *ofn;
        if(!fs::is_directory(output_dir))
        {
            spdlog::info(std::format("The directory provided does not exist. Creating directory {}", (std::string)output_dir));
            fs::create_directory(output_dir);
        }
    }

    std::size_t num_variations = program.get<std::size_t >("-N");
    std::size_t batch_size = program.get<std::size_t>("-b");
    double BBangle = program.get<double>("-q");
    double SCHangle = program.get<double>("-w");

    spdlog::info(std::format("Maximal BackBone angle set to {}. Maximal SideCHain angle set to {}.", BBangle, SCHangle));
    spdlog::info("Starting dataset generation.");

    createdataset(input_dir, output_dir, num_variations, batch_size, BBangle, SCHangle);


    spdlog::info("Dataset generation finished.");
    return 0;
}


