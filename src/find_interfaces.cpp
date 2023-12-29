#include<string>
#include<iostream>
#include<filesystem>
#include<array>
#include<vector>
#include<iterator>
#include<algorithm>
#include<thread>
#include<format>
#include<limits>
#include<cmath>
#include<chrono>
#include<argparse/argparse.hpp>
#include "molecules.h"
#include "io.h"
#include "fancy.h"
#include "threadpool.h"

//Verbose mode is currently not thread safe. I need to use mutexes or something...
using namespace std::chrono_literals;
namespace fs = std::filesystem;
bool verbose=false;
//bool force = false;

void createdataset(const std::string, const std::string, const unsigned int, const unsigned int, const bool);
void perturbRun(fs::path, fs::path, unsigned int, const bool);


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

    std::vector<fs::path> files = findInputFiles(input_dir);
    ProgressBar Pbar(files.size());

    unsigned i{0};
    for (auto& file:files)
    {
        if (!verbose) {
            i++;
            Pbar.update();
            std::string msg = std::to_string(static_cast<int>( i+1 )) + '/' +
                              std::to_string((int) (files.size()));
            Pbar.print(msg);
        }

        if (verbose){
            std::cout << "Opening " << file<< " for perturbation." << std::endl;
        }

        std::unique_ptr<RandomPerturbator> pert = std::make_unique<RandomPerturbator>(
                RandomPerturbator(file, verbose));

        if (verbose) {
            pert->getNumberOfResiduesPerChain(); //outputs how many residues there are in each chain.
        }
        if (verbose) std::cout << "Looking for interface residues." << std::endl;

        pert->findInterfaceResidues(12.0);

        if (verbose) std::cout << "Saving interface residues" << std::endl;
        std::string iface =  "_interfaces.json";
        std::string fname = file.stem();
        std::string fout = fname+iface;
        fs::path json_file{output_dir/fout};
        pert->saveInterfaceResidues(json_file);
        if (verbose) {
            pert->getInterfaceResidues();
        }
    }
    std::cout << "Dataset generation finished" << std::endl;
    return 0;
}

