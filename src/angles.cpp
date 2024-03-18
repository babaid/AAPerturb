#include "angles.h"
#include "spdlog/spdlog.h"

#include <random>
#include <vector>

//https://dunbrack.fccc.edu/bbdep2010/ConformationalAnalysis.php
double chi1() {
    auto logger = spdlog::get("main");
    // Seed the Mersenne Twister engine
    static std::random_device rd;
    static std::mt19937 mt(rd());

    // Define parameters for the normal distributions
    const std::vector<double> means = {120.0, 240.0, 360.0};
    const double stddev = 2.0;

    // Create the normal distributions
    std::vector<std::normal_distribution<double>> dists;
    for (const auto& mean : means) {
        dists.emplace_back(mean, stddev);
    }

    // Generate random numbers from the mixture
    std::uniform_int_distribution<> choicedist(0, 2);

    int choice = choicedist(mt);
    double ang = dists[choice](mt);
    logger->info(std::format("Chosen chi1 angle: {} degrees.", ang));
    return ang;
}

