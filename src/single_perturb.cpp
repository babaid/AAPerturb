//Main executable...
// Created by babaid on 07.09.23.
//
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
    program.add_argument("-N", "--num-variations")
            .scan<'d', std::size_t>()
            .default_value(std::size_t(1))
            .help("The number of variations of a single protein.");

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
    perturbRun(input_file, out, num_variations, verbose);
    //createdataset(input_dir, output_dir, num_variations, batch_size, verbose);
    std::cout << "Perturbation finished." << std::endl;
    return 0;
}





/*
 * Opens a PDB file and perturbes the interface amino acids in the protein a number of times.
 */
void perturbRun(fs::path input_filename, fs::path out,const unsigned int num_perturbations, const bool verbose) {
    std::size_t perturbcntr{number_of_files_in_directory(out)};
    if (perturbcntr<num_perturbations) {

        if (verbose){
            std::cout << "Opening " << input_filename << " for perturbation." << std::endl;
        }

        std::unique_ptr<RandomPerturbator> pert = std::make_unique<RandomPerturbator>(
                RandomPerturbator(input_filename, verbose));


        //auto path = out / "extracted";

        //if (!fs::is_directory(path)) fs::create_directory(path);
        //if (!fs::exists(path / "coords.tsv")) pert->saveCoords(path);

        if (verbose) {
            pert->getNumberOfResiduesPerChain(); //outputs how many residues there are in each chain.
        }
        if (verbose) std::cout << "Looking for interface residues." << std::endl;

        pert->findInterfaceResidues(12.0);

        if (verbose) {
            pert->getInterfaceResidues();
        }


        while (perturbcntr < num_perturbations) {

            std::string fname = std::to_string(perturbcntr) + ".pdb";
            fs::path out_path = out / fname;

            if (verbose) std::cout << "Choosing a random residue to perturb: ";

            std::pair<char, std::size_t> res = pert->chooseRandomResidue();

            if (verbose) std::cout << res.first << " : " << res.second << std::endl;
            Residue ref_residue = pert->getResidue(res.first, res.second);

            std::vector<std::string> comments;
            if (verbose) std::cout << "Perturbing the chosen residue";

            //Dont hate me but I get some random heap buffer overflow, so I will just deal with it later.
            double rmsd = std::numeric_limits<double>::infinity();
            try {
                pert->rotateResidueSidechainRandomly(res.first, res.second);
                pert->rotateResidueAroundBackboneRandomly(res.first, res.second);
                rmsd = pert->calculateRMSD(ref_residue);
            }
            catch (...) {
                if (verbose) std::cout << "Something was not right at " << input_filename << std::endl;
                continue;
            }
            if (verbose) std::cout << "Perturbation ended, per-residue RMSD: " << rmsd << std::endl;
            if (rmsd == 0.0) {
                perturbcntr--;
                continue;
            } else {
                perturbcntr++;
            }

            std::string comment1 = std::format("PERTURBATED RESIDUE: /{}:{}", res.first,
                                               std::to_string(res.second + 1));
            std::string comment2 = std::format("DISTMAT CHANGE INDEX : {} ",
                                               std::to_string(ref_residue.atoms.at(0).serial));
            std::string comment3 = std::format("RMSD: {}", rmsd);

            comments.push_back(comment1);
            comments.push_back(comment2);
            comments.push_back(comment3);

            if (verbose) std::cout << "Saving new PDB file at " << out_path << std::endl;
            //auto distmat_fname = std::to_string(perturbcntr) + ".tsv";
            // save update stuff.
            //pert->saveLocalDistMat(out / distmat_fname, res.first,
            //   res.second); //This either will make things fast or slow
            pert->saveToPDB(out_path, comments); //usual
            pert->setResidue(ref_residue); //hm.


        }
    }
}
