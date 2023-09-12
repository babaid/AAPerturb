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

void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein);
void perturbRun(fs::path filename, fs::path output_dir, unsigned int num_perturbations);

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("aaperturb");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output PDB which are perturbed randomly.");
    program.add_argument("-b", "--batch_size")
            .scan<'d', int>()
            .help("The number of threads to launch at dataset creation.");

    program.add_argument("-N", "--num-variations")
            .scan<'d', int>()
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
            std::cout << "The directory provided does not exist. Exiting... " << std::endl;
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
    std::size_t batch_size = 1;
    if(auto bsfn = program.present("-b"))
    {
        batch_size = std::stoi(*bsfn);
    }

    unsigned int num_variations = 10;
    if(auto numvar = program.present("-N"))
    {
        batch_size = std::stoi(*numvar);
    }



    createdataset(input_dir, output_dir, num_variations);

    return 0;
}





void perturbRun(fs::path filename, fs::path out, unsigned int num_perturbations) {

    std::map<char, std::vector<Residue *>> structure = parsePDB(filename);
    //std::map<char, std::vector<Residue*>> structure = parsePDB(filename);
    std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(structure, 9.0);
    for (unsigned int i = 0; i < num_perturbations; ++i) {
        std::string fname = std::to_string(i) + ".pdb";
        std::cout << "The output will be written to " << out / fname << std::endl;
        fs::path out_path = out / fname;

        std::pair<char, std::vector<std::size_t>> res = chooseRandomInterfaceResidue(structure,
                                                                                     interface_residue_indices);

        Residue *ref_residue = new Residue(*structure.at(res.first)[res.second[0]]);
        std::vector<std::string> comments;

        for (std::size_t &resid: res.second){
            rotateResidueSidechainRandomly(structure, res.first, resid);
            std::string comment = std::format("MUTATION: /{}:{}", res.first , std::to_string(resid+1));
            comments.push_back(comment);

        }



        saveToPDBWithComments(out_path, structure, comments);

        *structure.at(res.first)[res.second[0]] = *ref_residue;
        delete ref_residue;

    }
    for (const auto &chainEntry: structure) {
        for (Residue *residue: chainEntry.second) {
            delete residue;
        }

    }

}


void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein) {
    std::vector<fs::path> files = findInputFiles(inputdir);




    for (unsigned i=0; i<files.size(); ++i) {


            std::cout << "Opening " << files[i] << std::endl;
            fs::path filedir{files[i].filename()};
            filedir.replace_extension("");
            fs::path out = outputdir / filedir;
            fs::create_directory(out);
            perturbRun(files[i], out, num_variations_per_protein);



    }

}
