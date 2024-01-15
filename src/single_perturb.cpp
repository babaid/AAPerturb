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
            .help("Enable verbose mode. It is useful if you want to know if something is happening.")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output PDB which are perturbed randomly.");
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
    if(program["-v"] == true)
    {
        std::cout << "Running in verbose mode." << std::endl;
        verbose = true;
    }
    fs::path input_file, output_dir;
    if (auto ifn = program.present("-i"))
    {
        input_file = fs::absolute(*ifn);

//std::cout << "The directory provided does not exist. Exiting..." << std::endl;
  //          std::exit(1);

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
    double BBangle = program.get<double>("-q");
    double SCHangle = program.get<double>("-w");

    fs::path filedir{input_file.filename()};
    filedir.replace_extension("");
    fs::path out = output_dir / filedir;
    if (!fs::is_directory(out)) {
        fs::create_directory(out);
    }
    if (verbose) {

        std::cout << "Starting dataset generation." << std::endl;
        std::cout << "Output dir: " << out << std::endl;
    }
    perturbRun(input_file, out, num_variations, verbose, BBangle, SCHangle);
    //createdataset(input_dir, output_dir, num_variations, batch_size, verbose);
    std::cout << "Perturbation finished." << std::endl;
    return 0;
}




