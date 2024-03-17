#include <iostream>
#include <random>
#include <vector>
#include "angles.h"
//https://dunbrack.fccc.edu/bbdep2010/ConformationalAnalysis.php
double chi1() {

    // Seed the Mersenne Twister engine
    static std::random_device rd;
    static std::mt19937 mt(rd());

    // Define parameters for the normal distributions
    const std::vector<double> means = {120.0, 240.0, 360.0};
    const double stddev = 2.0;
    const std::vector<double> weights = {1/3, 1/3, 1/3};

    // Create the normal distributions
    std::vector<std::normal_distribution<double>> dists;
    for (const auto& mean : means) {
        dists.emplace_back(mean, stddev);
    }

    // Generate random numbers from the mixture
    double randomValue = mt() / static_cast<double>(mt.max()); // normalize to [0, 1]
    double cumulativeWeight = 0.0;

    for (size_t i = 0; i < means.size(); ++i) {
        cumulativeWeight += weights[i];
        if (randomValue < cumulativeWeight) {
            return dists[i](mt);
        }
    }
    // This line should not be reached under normal circumstances
    return dists.back()(mt);
}

