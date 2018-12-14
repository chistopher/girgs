
#include <girgs/Hyperbolic.h>

#include <random>
#include <fstream>

namespace girgs {


double calculateRadius(int n, double alpha, double T, int deg) {
    return 2 * log(n * 2 * alpha * alpha * (T == 0 ? 1 / M_PI : T / sin(M_PI * T)) /
                   (deg * (alpha - 0.5) * (alpha - 0.5)));
}

double hyperbolicDistance(double r1, double phi1, double r2, double phi2) {
    return acosh(std::max(1., cosh(r1 - r2) + (1. - cos(phi1 - phi2)) * sinh(r1) * sinh(r2)));
}

std::vector<double> sampleRadii(int n, double alpha, double R, int weightSeed) {
    auto result = std::vector<double>(n);
    auto gen = std::mt19937(weightSeed >= 0 ? weightSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i = 0; i < n; ++i) {
        auto p = dist(gen);
        while(p == 0) p = dist(gen);
        result[i] = acosh(p * (cosh(alpha * R) - 1) + 1) / alpha;
    }
    return result;
}

std::vector<double> sampleAngles(int n, int positionSeed) {
    auto result = std::vector<double>(n);
    auto gen = std::mt19937(positionSeed >= 0 ? positionSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i = 0; i < n; ++i)
        result[i] = dist(gen) * 2 * M_PI;
    return result;
}


} // namespace girgs
