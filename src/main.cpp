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
#include<argparse/argparse.hpp>
#include "../include/pdbparser.h"
#include "../include/geometry.h"
#include "../include/montecarlo.h"
#include "../include/io.h"

namespace fs = std::filesystem;

void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein, const unsigned int);
void perturbRun(fs::path filename, fs::path output_dir, unsigned int num_perturbations);

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("aaperturb");
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

    std::size_t num_variations = program.get<std::size_t >("-N");
    std::size_t batch_size = program.get<std::size_t>("-b");
    createdataset(input_dir, output_dir, num_variations, batch_size);
    return 0;
}





void perturbRun(fs::path filename, fs::path out, unsigned int num_perturbations) {

    std::map<char, std::vector<Residue *>>* structure = parsePDB(filename);
    //std::map<char, std::vector<Residue*>> structure = parsePDB(filename);
    std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(*structure, 9.0);
    for (unsigned int i = 0; i < num_perturbations; ++i) {
        std::string fname = std::to_string(i) + ".pdb";
        std::cout << "The output will be written to " << out / fname << std::endl;
        fs::path out_path = out / fname;

        if (!fs::exists(out_path)) {
            std::pair<char, std::vector<std::size_t>> res = chooseRandomInterfaceResidue(*structure,
                                                                                         interface_residue_indices);

            Residue *ref_residue = new Residue(*structure->at(res.first)[res.second[0]]);
            std::vector<std::string> comments;

            for (std::size_t &resid: res.second) {
                double rmsd = rotateResidueSidechainRandomly(*structure, res.first, resid);
                std::string comment1 = std::format("MUTATION: /{}:{}", res.first, std::to_string(ref_residue->resSeq ));
                std::string comment2 = std::format("RMSD: {}", rmsd);
                comments.push_back(comment1);
                comments.push_back(comment2);

            }


            saveToPDBWithComments(out_path, *structure, comments);

            *structure->at(res.first)[res.second[0]] = *ref_residue;
            delete ref_residue;
        }
        else
        {
            std::cout << "File already exists, skipping!" << std::endl;
        }

    }

    for (const auto &chainEntry: *structure)
    {
        for (Residue *residue: chainEntry.second)
        {
            delete residue;
        }

    }
    delete structure;

}


void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein, const unsigned int batch_size ) {
    std::vector<fs::path> files = findInputFiles(inputdir);


    for (unsigned int i{0}; i<files.size(); i++) {
            std::vector<std::thread> ThreadVector;
            while((i%batch_size != 0 || i==0) && i<files.size()) {
                std::cout << "Opening " << files[i] << std::endl;
                fs::path filedir{files[i].filename()};
                filedir.replace_extension("");
                fs::path out = outputdir / filedir;
                fs::create_directory(out);
                ThreadVector.emplace_back(perturbRun, files[i], out, num_variations_per_protein);
                ++i;
                //perturbRun(files[i], out, num_variations_per_protein)
            }
            for(auto& t:ThreadVector) t.join();
            ThreadVector.clear();

    }

}
