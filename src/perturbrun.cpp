#include<string>
#include<iostream>
#include<filesystem>
#include<array>
#include<vector>
#include<algorithm>
#include<thread>
#include<format>
#include<limits>
#include<argparse/argparse.hpp>
#include "molecules.h"
#include "io.h"
#include "fancy.h"
#include "threadpool.h"
#include "perturbrun.h"


using namespace std::chrono_literals;
namespace fs = std::filesystem;





/*
 * Opens a PDB file and perturbes the interface amino acids in the protein a number of times.
 */
void perturbRun(fs::path input_filename, fs::path out,const unsigned int num_perturbations, const bool verbose, double BBangle, double SCHangle) {
    std::size_t perturbcntr{number_of_files_in_directory(out)};
    if (perturbcntr<num_perturbations) {

        if (verbose){
            std::cout << "Opening " << input_filename << " for perturbation." << std::endl;
        }

        std::unique_ptr<RandomPerturbator> pert = std::make_unique<RandomPerturbator>(
                RandomPerturbator(input_filename, verbose));

        pert->setMaxRotAngleBB(BBangle);
        pert->setMaxRotAngleSCH(SCHangle);

        //auto path = out / "extracted";

        //if (!fs::is_directory(path)) fs::create_directory(path);
        //if (!fs::exists(path / "coords.tsv")) pert->saveCoords(path);

        if (verbose) {
            pert->getNumberOfResiduesPerChain(); //outputs how many residues there are in each chain.
        }
        if (verbose) std::cout << "Looking for interface residues." << std::endl;

        pert->findInterfaceResidues(12.0);

        if (verbose) std::cout << "Saving interface residues" << std::endl;
        fs::path json_file{out /"interfaces.json"};
        pert->saveInterfaceResidues(json_file);
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
        pert.reset();
    }

}


void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein, const unsigned int batch_size, const bool verbose, double BBangle, double SCHangle) {

    //Batched threadpool. We wait calmly for each thread to finish, when they finished, we empty the queue of tasks, and start the next iteration.
    // I strictly want to enqueue just as many tasks as the batch size allows, avoiding any type of overflows
    ThreadPool pool(batch_size); // Thread pool UwU
    std::vector<std::future<void>> futures;
    std::vector<fs::path> files = findInputFiles(inputdir);
    ProgressBar Pbar(files.size());


    if(!verbose){Pbar.print("0/0");}
    for (unsigned int batch_start{0}; batch_start < files.size();batch_start+=batch_size) {
        std::this_thread::sleep_for(0.5s);
        if (!verbose) {
            Pbar.update();
            std::string msg = std::to_string(static_cast<int>( batch_start / batch_size + 1 )) + '/' +
                              std::to_string((int) (files.size() / batch_size));
            Pbar.print(msg);
        }
        else {
            std::cout << "Working on batch " << static_cast<int>(batch_start / batch_size + 1 )
                      << "/" << (int) (files.size() / batch_size)
                      << "." << std::endl;
        }

        unsigned int batch_end = std::min(batch_start + batch_size, static_cast<unsigned int>(files.size()));
        for(unsigned int i{batch_start}; i<batch_end;++i){
            //filesystem stuff
            fs::path filedir{files[i].filename()};
            filedir.replace_extension("");
            fs::path out = outputdir / filedir;
            if (fs::is_directory(out)) {
                if (number_of_files_in_directory(out) == num_variations_per_protein) {}//Pbar.update();}
            } else fs::create_directory(out);


            // filesystem stuff done
            std::future<void> result = pool.enqueue(perturbRun, files[i], out, num_variations_per_protein, verbose, BBangle, SCHangle);
            futures.emplace_back(std::move(result));

        }


        //std::this_thread::sleep_for(std::chrono::seconds((int)batch_size)); //longest operation takes about a second

        futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
            return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        }), futures.end());

        //Not sure if this is good practice but it avoids enqueuing too much stuff
        if(futures.size()>batch_size*10){
            while(futures.size()!=0) {
                std::this_thread::sleep_for(std::chrono::seconds(batch_size)); //wait for five seconds so tasks can finish
                futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
                    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
                }), futures.end()); // delete finished tasks
            }
        } //wait for some threads to finish, so we dont overload
    }
}
