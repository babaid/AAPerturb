#include "molecules.h"
#include "IOutils.h"
#include "fancy.h"
#include "threadpool.h"
#include "perturbrun.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include<string>
#include<filesystem>
#include<array>
#include<vector>
#include<algorithm>
#include<thread>
#include<format>
#include<limits>





using namespace std::chrono_literals;
namespace fs = std::filesystem;


/*
 * Opens a PDB file and perturbes the interface amino acids in the protein a number of times.
 */
void perturbRun(fs::path& input_filename, fs::path& out, size_t num_perturbations, double BBangle, double SCHangle) {
    auto console = spdlog::get("main");

    std::size_t  cyclecntr{0}, perturbcntr{0};

    if (perturbcntr<num_perturbations) {
        console->info(std::format("Opening {} for perturbation procedure.", input_filename.filename().string()));

        std::unique_ptr<RandomPerturbator> pert = std::make_unique<RandomPerturbator>(RandomPerturbator(input_filename, BBangle, SCHangle));

        pert->getNumberOfResiduesPerChain(); //outputs how many residues there are in each chain.
        pert->findInterfaceResidues(12.0);
        fs::path json_file{out / "interfaces.json"};
        pert->saveInterfaceResidues(json_file);
        pert->printInterfaceResidues();
        auto interface_residues = pert->getInterfaceResidues();

        while (cyclecntr < num_perturbations) {
            for (const auto &entry: interface_residues) {
                char chain = entry.first;
                for (size_t resid : entry.second) {
                    std::string fname = std::to_string(perturbcntr) + ".pdb";
                    fs::path out_path = out / fname;

                    console->info(std::format("Perturbing residue {}/{}", entry.second.size(), perturbcntr));
                    console->info(std::format("{} : {}", chain, resid));


                    Residue ref_residue = pert->getResidue(chain, resid);

                    console->info("Started perturbation calculations.");

                    std::vector<std::string> comments;
                    double rmsd = std::numeric_limits<double>::infinity();

                    try {
                        pert->rotateResidueSidechainRandomly(chain, resid);
                        pert->rotateResidueAroundBackboneRandomly(chain, resid);
                        rmsd = pert->calculateRMSD(ref_residue);
                    }
                    catch (...) {
                        console->critical(std::format("Something was not right at {}. Skipping.", input_filename.string()));
                        continue;
                    }

                    console->info(std::format("Perturbation done. Per-Resiude RMSD: {}", rmsd));

                    if (rmsd == 0.0) {
                        perturbcntr--;
                        continue;
                    } else {
                        perturbcntr++;
                    }

                    std::string comment1 = std::format("PERTURBATED RESIDUE: /{}:{}", chain,
                                                       std::to_string(resid + 1));
                    std::string comment2 = std::format("DISTMAT CHANGE INDEX : {} ",
                                                       std::to_string(ref_residue.atoms.at(0).serial));
                    std::string comment3 = std::format("RMSD: {}", rmsd);

                    comments.push_back(comment1);
                    comments.push_back(comment2);
                    comments.push_back(comment3);
                    pert->saveToPDB(out_path, comments); //usual
                    pert->setResidue(ref_residue); //hm.
                }
            }
        cyclecntr++;
        }
    }

}


void createdataset(const std::string& inputdir, const std::string& outputdir, size_t num_variations_per_protein, size_t batch_size, double BBangle, double SCHangle) {
    auto console = spdlog::get("main");
    //Batched threadpool. We wait calmly for each thread to finish, when they finished, we empty the queue of tasks, and start the next iteration.
    // I strictly want to enqueue just as many tasks as the batch size allows, avoiding any type of overflows
    ThreadPool pool(batch_size); // Thread pool UwU
    std::vector<std::future<void>> futures;
    std::vector<fs::path> files = findInputFiles(inputdir);

    ProgressBar Pbar((int)files.size());

    if(console->level() == spdlog::level::off) Pbar.print("0/0");

    for (size_t batch_start{0}; batch_start < files.size();batch_start+=batch_size) {
        std::this_thread::sleep_for(0.1s);
        ///
        if(console->level() == spdlog::level::off)
        {
            Pbar.update();
            std::string msg = std::to_string(static_cast<int>( batch_start / batch_size + 1 )) + '/' +
                              std::to_string((int) (files.size() / batch_size));

            Pbar.print(msg);
        }
        else console->info(std::format("Working on batch {}/{}.", static_cast<int>(batch_start / batch_size + 1 ), (int) (files.size() / batch_size)));
        ////

        size_t batch_end = std::min(unsigned(batch_start + batch_size), static_cast<unsigned int>(files.size()));
        for(size_t i{batch_start}; i<batch_end;++i){

            //filesystem stuff
            fs::path filedir{files[i].filename()};
            filedir.replace_extension("");
            fs::path out = outputdir / filedir;
            if (!fs::is_directory(out)) fs::create_directory(out);
            // filesystem stuff done

            std::future<void> result = pool.enqueue(perturbRun, files[i], out, num_variations_per_protein, BBangle, SCHangle);
            futures.emplace_back(std::move(result));

        }

        //std::this_thread::sleep_for(std::chrono::seconds((int)batch_size)); //longest operation takes about a second
        futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
            return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        }), futures.end());

        //Not sure if this is good practice but it avoids enqueuing too much stuff
        if(futures.size()>batch_size*10){
            while(!futures.empty()) {
                std::this_thread::sleep_for(std::chrono::seconds(batch_size)); //wait for five seconds so tasks can finish
                futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
                    return f.wait_for(std::chrono::milliseconds(100)) == std::future_status::ready;
                }), futures.end()); // delete finished tasks
            }
        } //wait for some threads to finish
    }
}
