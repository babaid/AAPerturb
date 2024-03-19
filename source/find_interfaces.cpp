#include "molecules.h"
#include "IOutils.h"
#include "fancy.h"
#include "threadpool.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include<argparse/argparse.hpp>

#include<string>
#include<iostream>
#include<filesystem>
#include<vector>
#include<iterator>
#include<cmath>


//Verbose mode is currently not thread safe. I need to use mutexes or something...
using namespace std::chrono_literals;
namespace fs = std::filesystem;


int main(int argc, char *argv[]) {
    auto logger = spdlog::stdout_logger_mt("main");
    argparse::ArgumentParser program("find_interfaces", "1.2.0");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output interface json files.");

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
            logger->info(std::format("The directory provided does not exist. Creating directory ", output_dir.string()));
            fs::create_directory(output_dir);
        }
    }

    std::vector<fs::path> files = findInputFiles(input_dir);
    ProgressBar Pbar(files.size());

    unsigned i{0};
    for (auto& file:files)
    {
        if (logger->level() == spdlog::level::off) {
            i++;
            Pbar.update();
            std::string msg = std::to_string(static_cast<int>( i+1 )) + '/' +
                              std::to_string((int) (files.size()));
            Pbar.print(msg);
        }

        logger->info(std::format("Opening {} for perturbation.", file.string()));

        std::unique_ptr<RandomPerturbator> pert = std::make_unique<RandomPerturbator>(
                RandomPerturbator(file));


        pert->getNumberOfResiduesPerChain(); //outputs how many residues there are in each chain.

        logger->info("Looking for interface residues.");

        pert->findInterfaceResidues(12.0);

        logger->info("Saving interface residues.");

        std::string iface =  "_interfaces.json";
        std::string fname = file.stem();
        std::string fout = fname+iface;
        fs::path json_file{output_dir/fout};
        pert->saveInterfaceResidues(json_file);
        pert->getInterfaceResidues();

    }
    logger->info("Interface lookup finished.");
    return 0;
}

