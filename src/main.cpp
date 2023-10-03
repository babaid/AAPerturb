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
#include "pdbparser.h"
#include "geometry.h"
#include "montecarlo.h"
#include "io.h"
#include "fancy.h"
#include "threadpool.h"

//Verbose mode is currently not thread safe. I need to use mutexes or something...
using namespace std::chrono_literals;
namespace fs = std::filesystem;       
bool verbose=false;
bool force = false;

void createdataset(const std::string, const std::string, const unsigned int, const unsigned int, const bool, const bool);
void perturbRun(fs::path, fs::path, unsigned int, const bool);
std::size_t number_of_files_in_directory(fs::path path);
int main(int argc, char *argv[]) {

    argparse::ArgumentParser program("aaperturb");
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
    program.add_argument("-b", "--batch-size")
            .scan<'d', std::size_t>()
            .default_value(std::size_t(1))
            .help("The number of batches for multiprocessing");
    program.add_argument("-N", "--num-variations")
            .scan<'d', std::size_t>()
            .default_value(std::size_t(1))
            .help("The number of variations of a single protein.");
    program.add_argument("-f", "--force")
            .help("Force recreation of already existent files in the output directory. Treat with care.")
            .default_value(false)
            .implicit_value(true);

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

    if(program["-f"] == true)
    {
        if (verbose) std::cout << "You forced me to use force." << std::endl;
        force=true;
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


    //auto pdb = parsePDB("../_out/1qse.pdb", verbose=true);
    /*for (auto& chain:*pdb)
    {
        std::cout << chain.first << ":" << chain.second.size() << std::endl;
        if (chain.first == 'D') {
            for (auto &el: chain.second) std::cout << el.resSeq << std::endl;
        }
    }*/

    std::cout<< "Starting"<< std::endl;
    createdataset(input_dir, output_dir, num_variations, batch_size, force, verbose);
    std::cout << "Dataset creation finished" << std::endl;
    return 0;
}





/*
 * Opens a PDB file and perturbes the interface amino acids in the protein a number of times.
 */
void perturbRun(fs::path filename, fs::path out,const unsigned int num_perturbations, const bool verbose) {

    if (verbose){
        std::cout << "Opening " << filename << " for perturbation." << std::endl;
    }
    std::unique_ptr<std::map<char, std::vector<Residue>>>  structure = parsePDB(filename);

    if (verbose)
    {
        auto print_chain_n = [](auto const& elem){std::cout << elem.first << ": " << elem.second.size() << ", ";};
        std::for_each(structure->begin(), structure->end(), print_chain_n);
    }
    if (verbose) std::cout << "Looking for interface residues." << std::endl;

    std::map<char, std::vector<int>> interface_residue_indices = findInterfaceResidues(structure, 12.0);

    if (verbose){
        std::cout << "Found following number of interface residues in the chains: ";
        auto print_chain_n = [](auto const& elem){std::cout << elem.first << ": " << elem.second.size() << ", ";};
        std::for_each(interface_residue_indices.begin(), interface_residue_indices.end(), print_chain_n);
        std::cout<< std::endl;
    }

    std::size_t perturbcntr{number_of_files_in_directory(out)};

    while(perturbcntr<num_perturbations) {

        std::string fname = std::to_string(perturbcntr) + ".pdb";
        fs::path out_path = out / fname;

        if(verbose) std::cout <<  "Choosing a random residue to perturb: ";

        std::pair<char , std::size_t>  res = chooseRandomResidue(interface_residue_indices);

        if(verbose) std::cout << res.first << " : " << res.second << std::endl;

        Residue ref_residue(structure->at(res.first).at(res.second));
        std::vector<std::string> comments;



        if(verbose) std::cout << "Perturbing the chosen residue";

        //Dont hate me but I get some random heap buffer overflow, so I will just deal with it later.
        double rmsd = std::numeric_limits<double>::infinity();
        try {
            rmsd = rotateResidueSidechainRandomly(structure, res.first, res.second, verbose);
        }
        catch (...)
        {
            if(verbose)std::cout << "Something was not right at "  << filename <<std::endl;
            continue;
        }

        if(verbose) std::cout << "Perturbation ended, per-residue RMSD: " << rmsd << std::endl;
        if (rmsd == 0.0){
            perturbcntr--; continue;
        } else {
            perturbcntr++;
        }

            std::string comment1 = std::format("MUTATION: /{}:{}", res.first, std::to_string(ref_residue.resSeq));
            std::string comment2 = std::format("RMSD: {}", rmsd);

            comments.push_back(comment1);
            comments.push_back(comment2);

            if (verbose) std::cout << "Saving new PDB file at " << out_path << std::endl;

            saveToPDBWithComments(out_path, structure, comments);
            structure->at(res.first).at(res.second) = ref_residue;

    }
}


void createdataset(const std::string inputdir, const std::string outputdir, const unsigned int num_variations_per_protein, const unsigned int batch_size,  const bool force, const bool verbose) {

    //Batched threadpool. We wait calmly for each thread to finish, when they finished, we empty the queue of tasks, and start the next iteration.
    // I strictly want to enqueue just as many tasks as the batch size allows, avoiding any type of overflows
    ThreadPool pool(batch_size); // Thread pool UwU
    std::vector<std::future<void>> futures;
    std::vector<fs::path> files = findInputFiles(inputdir);
    ProgressBar Pbar(files.size());
    Pbar.print("0/0");
    for (unsigned int batch_start{0}; batch_start < files.size();batch_start+=batch_size) {
        int batch_end = std::min(batch_start + batch_size, static_cast<unsigned int>(files.size()));
        for(unsigned int i{batch_start}; i<batch_end;++i){
            //filesystem stuff
            fs::path filedir{files[i].filename()};
            filedir.replace_extension("");
            fs::path out = outputdir / filedir;
            if (fs::is_directory(out)) {
                if (!verbose) {
                    std::string msg = std::to_string(static_cast<int>(i / batch_size + 1)) + '/' +
                                      std::to_string((int) (files.size() / batch_size + 1));
                    Pbar.print(msg);
                }
                if (number_of_files_in_directory(out) == num_variations_per_protein) {}//Pbar.update();}
            } else fs::create_directory(out);


            // filesystem stuff done
            std::future<void> result = pool.enqueue(perturbRun, files[i], out, num_variations_per_protein, verbose);
            futures.emplace_back(std::move(result));

            if (!verbose) {

                Pbar.update();
                std::string msg = std::to_string(static_cast<int>(i / batch_size + 1 )) + '/' +
                                  std::to_string((int) (files.size() / batch_size + 1));
                Pbar.print(msg);
            } else  std::cout << "Batch " << static_cast<int>(i / batch_size + 1 )
                              << "/" << (int) (files.size() / batch_size + 1)
                              << " is done." << std::endl;
        }

        std::this_thread::sleep_for(std::chrono::seconds((int)batch_size/2)); //longest operation takes about a second

        futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
            return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        }), futures.end());

        //Not sure if this is good practice but it avoids enqueuing too much stuff
        if(futures.size()>batch_size*2){
            while(futures.size()!=0) {
                std::this_thread::sleep_for(std::chrono::seconds(batch_size)); //wait for five seconds so tasks can finish
                futures.erase(std::remove_if(futures.begin(), futures.end(), [](const std::future<void> &f) {
                    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
                }), futures.end()); // delete finished tasks
            }
        } //wait for some threads to finish, so we dont overload
    }
}


std::size_t number_of_files_in_directory(fs::path path)
{
    using std::filesystem::directory_iterator;
    using fp = bool (*)( const std::filesystem::path&);
    return std::count_if(directory_iterator(path), directory_iterator{}, (fp)std::filesystem::is_regular_file);
}