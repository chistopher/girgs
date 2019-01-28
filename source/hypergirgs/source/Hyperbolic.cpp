
#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/HyperbolicTree.h>

#include <random>
#include <fstream>
#include <cmath>
#include <mutex>

#include <omp.h>


namespace hypergirgs {


double calculateRadius(int n, double alpha, double T, int deg) {
    return 2 * log(n * 2 * alpha * alpha * (T == 0 ? 1 / PI : T / sin(PI * T)) /
                   (deg * (alpha - 0.5) * (alpha - 0.5)));
}

double hyperbolicDistance(double r1, double phi1, double r2, double phi2) {
    return acosh(std::max(1., cosh(r1 - r2) + (1. - cos(phi1 - phi2)) * sinh(r1) * sinh(r2)));
}

std::vector<double> sampleRadii(int n, double alpha, double R, int seed, bool parallel) {
    std::vector<double> result(n);

    auto threads = parallel ? omp_get_max_threads() : 1;
    auto gens = std::vector<hypergirgs::default_random_engine>(threads);
    auto dists = std::vector<std::uniform_real_distribution<>>(threads);
    for (int thread = 0; thread < threads; thread++) {
        gens[thread].seed(seed >= 0 ? seed + thread : std::random_device()());
    }

    const auto invalpha = 1.0 / alpha;
    const auto factor = std::cosh(alpha * R) - 1.0;

    #pragma omp parallel num_threads(threads)
    {
        auto& gen = gens[omp_get_thread_num()];
        auto& dist = dists[omp_get_thread_num()];
        #pragma omp for
        for (int i = 0; i < n; ++i) {
            auto p = dist(gen);
            while (p == 0) p = dist(gen);
            result[i] = acosh(p * factor + 1.0) * invalpha;
        }
    }

    return result;
}

std::vector<double> sampleAngles(int n, int seed, bool parallel) {
    std::vector<double> result(n);
    
    auto threads = parallel ? omp_get_max_threads() : 1;
    auto gens = std::vector<hypergirgs::default_random_engine>(threads);
    auto dists = std::vector<std::uniform_real_distribution<>>();
    for (int thread = 0; thread < threads; thread++) {
        gens[thread].seed(seed >= 0 ? seed + thread : std::random_device()());
        dists.emplace_back(0.0, std::nextafter(2 * PI, 0.0));
    }

    #pragma omp parallel num_threads(threads)
    {
        auto& gen = gens[omp_get_thread_num()];
        auto& dist = dists[omp_get_thread_num()];
        #pragma omp for
        for (int i = 0; i < n; ++i)
            result[i] = dist(gen);
    }

    return result;
}

std::vector<std::pair<int, int> > generateEdges(std::vector<double>& radii, std::vector<double>& angles, double T, double R, int seed) {
    std::vector<std::pair<int,int>> graph;
    std::mutex m;

    auto addEdge = [&graph, &m] (int u, int v, int tid) {
        std::lock_guard<std::mutex> lock(m);
        graph.emplace_back(u,v);
    };

    auto generator = hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge);
    generator.generate(seed);

    return graph;
}

} // namespace hypergirgs
