//Main executable...
// Created by babaid on 07.09.23.
//

#include<iostream>
#include<filesystem>
#include<iterator>
#include<thread>
#include<cmath>
#include<argparse/argparse.hpp>
#include "perturbrun.h"

//Verbose mode is currently not thread safe. I need to use mutexes or something...
using namespace std::chrono_literals;
namespace fs = std::filesystem;
bool verbose=false;
//bool force = false;




int main(int argc, char *argv[]) {

    argparse::ArgumentParser program("aaperturb", "1.2.0");
    program.add_argument("-v", "--verbose")
            .help("Enable verbose mode. You should use this only with a batch size of 1, otherwise weird stuff could happen.")
            .default_value(false)
            .implicit_value(true);
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
    //program.add_argument("-f", "--force")
    //        .help("Force recreation of already existent files in the output directory. Treat with care.")
    //        .default_value(false)
    //        .implicit_value(true);

    try {
        program.parse_args(argc, argv);
        }
    catch (const std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }
    if(program["-v"] == true)
    {
        std::cout << "Running in verbose mode." << std::endl;
        verbose = true;
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


    std::size_t num_variations = program.get<std::size_t >("-N");
    std::size_t batch_size = program.get<std::size_t>("-b");
    double BBangle = program.get<double>("-q");
    double SCHangle = program.get<double>("-w");

    if (verbose) std::cout << "Maximal BB angle set to: " << BBangle<< std::endl<< "Maximal SCH angle set to: " << SCHangle << std::endl;
    std::cout<< "Starting dataset generation."<< std::endl;
    createdataset(input_dir, output_dir, num_variations, batch_size, verbose, BBangle, SCHangle);
    std::cout << "Dataset generation finished" << std::endl;
    return 0;
}


