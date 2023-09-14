//
// Created by babaid on 14.09.23.
//
#include<string>
#include<iostream>
#include<filesystem>
#include<vector>
#include<iterator>
#include<algorithm>
#include<thread>
#include<argparse/argparse.hpp>
#include "../include/pdbparser.h"
#include "../include/io.h"

class pBar {
public:
    void update(double newProgress) {
        currentProgress += newProgress;
        amountOfFiller = (int)((currentProgress / neededProgress)*(double)pBarLength);
    }
    void print() {
        currUpdateVal %= pBarUpdater.length();
        std::cout << "\r" //Bring cursor to start of line
             << firstPartOfpBar; //Print out first part of pBar
        for (int a = 0; a < amountOfFiller; a++) { //Print out current progress
            std::cout << pBarFiller;
        }
        std::cout << pBarUpdater[currUpdateVal];
        for (int b = 0; b < pBarLength - amountOfFiller; b++) { //Print out spaces
            std::cout << " ";
        }
        std::cout << lastPartOfpBar //Print out last part of progress bar
             << " (" << (int)(100*(currentProgress/neededProgress)) << "%)" //This just prints out the percent
             << std::flush;
        currUpdateVal += 1;
    }
    std::string firstPartOfpBar = "[", //Change these at will (that is why I made them public)
    lastPartOfpBar = "]",
            pBarFiller = "|",
            pBarUpdater = "/-\\|";
private:
    int amountOfFiller,
            pBarLength = 50, //I would recommend NOT changing this
    currUpdateVal = 0; //Do not change
    double currentProgress = 0, //Do not change
    neededProgress = 100; //I would recommend NOT changing this
};


int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("pdbcleaner");
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

    pBar bar;
    auto files = findInputFiles(input_dir);
    for (std::size_t i{0}; i<files.size();++i) {
        bar.update(100/files.size());
        if (!fs::exists(output_dir/files[i].filename())) {
            auto clean_structure = parsePDBToBeCleaned(files[i]);
            std::vector<std::string> comments;
            comments.push_back("This file was previously reindexed and the waters and the hydrogens were removed");
            saveToPDBWithComments(output_dir / files[i].filename(), clean_structure, comments);
        }
        bar.print();
    }
    return 0;
}



