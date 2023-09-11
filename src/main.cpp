//Main executable...
// Created by babaid on 07.09.23.
//
#include<string>
#include<iostream>
#include<filesystem>
#include<array>
#include<vector>
#include<random>
#include<iterator>
#include<algorithm>
#include<thread>
#include<functional>
#include<argparse/argparse.hpp>
#include "../include/pdbparser.h"
#include "../include/geometry.h"
#include "../include/montecarlo.h"
#include "../include/molecules.h"
#include "../include/io.h"
namespace fs = std::filesystem;

void createdataset(const std::string inputdir, const std::string outputdir, const std::size_t batch_size, const unsigned int num_variations_per_protein);
void perturbRun(fs::path filename, fs::path output_dir, unsigned int num_perturbations);

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("aasp");
    program.add_argument("-i", "--input-dir")
            .required()
            .help("The directory containing the input PDB files which we want to perturb randomly.");
    program.add_argument("-o", "--output-dir")
            .required()
            .help("The directory containing the output PDB which are perturbed randomly.");
    program.add_argument("-b", "--batch_size")
            .scan<'d', int>()
            .help("The number of threads to launch at dataset creation.");

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
            std::exit(1);
        }


    }
    std::size_t batch_size = 1;
    if(auto bsfn = program.present("-b"))
    {
        batch_size = std::stoi(*bsfn);
    }

    std::map<char, std::vector<Residue*>> structure = parsePDB("_data/1nwo.pdb");
    //std::map<char, std::vector<int>> ifr = findInterfaceResidues(structure);
    for(Residue* res: structure['A'])

    {
        res->atoms[0].coords[0] += 1;
        std::cout << res->resName <<" "<< res->atoms[0].name << " LINE: "<< res->atoms[0].serial<< std::endl;
    }
    saveToPDB("test.pdb", structure);

    std::cout << std::endl;
    std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(structure);
    for (auto el: interface_residue_indices['A'])
    {
        std::cout << el << std::endl;
    }

    createdataset(input_dir, output_dir, batch_size, 5);

}





void perturbRun(fs::path filename, fs::path out, unsigned int num_perturbations)
{

    std::map<char, std::vector<Residue*>>* structure = new std::map<char, std::vector<Residue*>>(parsePDB(filename));
    //std::map<char, std::vector<Residue*>> structure = parsePDB(filename);
    std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(*structure);

    for(unsigned int i =0; i<num_perturbations;++i)
    {
        std::string fname = std::to_string(i) + ".pdb";
        std::cout << "The output will be written to " << out/fname << std::endl;

        std::pair<char, std::vector<std::size_t>> res = chooseRandomInterfaceResidue(*structure, interface_residue_indices);

        Residue* ref_residue = new Residue(*structure->at(res.first)[res.second[0]]);

        for(std::size_t& resid:res.second) rotateResidueSidechainRandomly(*structure, res.first, resid);

        saveToPDB(out/fname, *structure);

        *structure->at(res.first)[res.second[0]] = *ref_residue;

    }
}

void createdataset(const std::string inputdir, const std::string outputdir, const std::size_t batch_size, const unsigned int num_variations_per_protein) {
    std::vector<fs::path> files = createFileBatches(inputdir, batch_size);


    for (auto &filename: files) {
        //std::vector<std::thread> ThreadVector;
        std::cout  << "Opening " << filename << std::endl;
        fs::path filedir{filename.filename()};
        filedir.replace_extension("");
        fs::path out = outputdir / filedir;
        fs::create_directory(out);
        perturbRun(filename, out, num_variations_per_protein);
        //T//hreadVector.emplace_back(std::thread([&](){ perturbRun(filename, out, num_variations_per_protein);}));

        //std::for_each(ThreadVector.begin(), ThreadVector.end(), [](std::thread &t){t.join();});
    }

}
