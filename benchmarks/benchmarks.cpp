//
// Created by babaid on 01.10.23.
//

#include<argparse/argparse.hpp>
#include "pdbparser.h"
#include "io.h"
#include "geometry.h"
#include "montecarlo.h"
#include<filesystem>
#include<chrono>
#include<string>
#include<vector>
#include "benchmarks.h"


namespace fs = std::filesystem;

void benchmark_pdbparser(fs::path& input_path,fs::path& output_path, unsigned cnt)
{

    auto start = std::chrono::high_resolution_clock::now();
    auto files = findInputFiles(input_path);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time spent on finding all the input files: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start) << std::endl;
    std::chrono::milliseconds read_duration{0};
    std::chrono::milliseconds write_duration{0};
    std::chrono::milliseconds iffind_duration{0};
    std::chrono::milliseconds reschoose_duration{0};
    std::chrono::milliseconds resrot_duration{0};
    for (unsigned i{0}; i<cnt; i++){
        for (auto& file:files) {
            start = std::chrono::high_resolution_clock::now();
            auto structure = parsePDB(file);
            end = std::chrono::high_resolution_clock::now();
            read_duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);

            start = std::chrono::high_resolution_clock::now();
            auto interfaces = findInterfaceResidues(structure, 9.0);
            end = std::chrono::high_resolution_clock::now();
            iffind_duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);

            start = std::chrono::high_resolution_clock::now();
            auto ifres = chooseRandomResidue(interfaces);
            end = std::chrono::high_resolution_clock::now();
            reschoose_duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);

            start = std::chrono::high_resolution_clock::now();

            double rmsd = rotateResidueSidechainRandomly(structure, ifres.first, ifres.second, false);

            end = std::chrono::high_resolution_clock::now();
            resrot_duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);

            std::string filename = std::to_string(i)+ "_" + (std::string)file.filename();

            std::string c1 = "RMSD: " + std::to_string(rmsd);
            std::string c2 = "Rotated Residue: " + std::to_string(ifres.first) + '/' + std::to_string(ifres.second);
            std::vector<std::string> comments{c1, c2};

            start = std::chrono::high_resolution_clock::now();
            saveToPDBWithComments(output_path/filename, structure, comments);
            end = std::chrono::high_resolution_clock::now();
            write_duration+=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        }
        read_duration/=files.size();
        write_duration/=files.size();
        iffind_duration/=files.size();
        reschoose_duration/=files.size();
        resrot_duration/=files.size();
    }
    read_duration/=(double)cnt;
    write_duration/=(double)cnt;
    iffind_duration/=(double)cnt;
    reschoose_duration/=(double)cnt;
    resrot_duration/=(double)cnt;


    std::cout << "Overall avg. time for parsing a PDB file: " << read_duration << std::endl;
    std::cout << "Overall avg. time for saving a PDB file: " << write_duration << std::endl;
    std::cout << "Overall avg. time for finding the interfaces: " << iffind_duration << std::endl;
    std::cout << "Overall avg. time for choosing a residue: " << reschoose_duration << std::endl;
    std::cout << "Overall avg. time for rotating a residue: " << resrot_duration << std::endl;



    std::cout << "Approximated time for a dataset of 3e3 raw structures, with 1e5 varitions each: "
                << std::chrono::duration_cast<std::chrono::minutes>((read_duration+write_duration+iffind_duration)*3000+ (reschoose_duration+resrot_duration)*10000)
                        << std::endl << "This does not include possible 0 RMSD perturbations.";

}
