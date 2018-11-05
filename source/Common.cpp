
#include <Common.h>

#include <random>
#include <cmath>


std::vector<double> generatePowerLawWeights(unsigned int n, double lower, double upper, double beta, int seed = -1) {
    auto gen = std::mt19937(seed >= 0 ?  seed : std::random_device()());
    std::uniform_real_distribution<> dist(0,1);
    auto weights = std::vector<double>(n);
    for(auto i=0; i<n; ++i)
        weights[i] = std::pow(
                (std::pow(upper,beta+1)-std::pow(lower,beta+1))*dist(gen) + std::pow(lower,beta+1),
                1.0/(beta+1) );
    return weights;
}


std::vector<double> generateWeights(unsigned int n, double beta, int seed) {
    return generatePowerLawWeights(n, 1.0, n, beta, seed);
}

