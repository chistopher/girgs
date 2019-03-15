
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

    const auto threads = parallel ? omp_get_max_threads() : 1;

    const auto invalpha = 1.0 / alpha;
    #pragma omp parallel num_threads(threads)
    {
        auto gen = hypergirgs::default_random_engine(
            seed >= 0 ? seed + omp_get_thread_num() : std::random_device()());
        auto dist = std::uniform_real_distribution<>(std::nextafter(1.0, 2.0), std::cosh(alpha * R));
        
        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            result[i] = acosh(dist(gen)) * invalpha;
        }
    }

    return result;
}

std::vector<double> sampleAngles(int n, int seed, bool parallel) {
    std::vector<double> result(n);
    
    const auto threads = parallel ? omp_get_max_threads() : 1;

    #pragma omp parallel num_threads(threads)
    {
        auto gen = hypergirgs::default_random_engine(
            seed >= 0 ? seed + omp_get_thread_num() : std::random_device()());
        auto dist = std::uniform_real_distribution<>(0, 2 * PI);

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i)
            result[i] = dist(gen);
    }

    return result;
}

std::vector<std::pair<int, int> > generateEdges(std::vector<double>& radii, std::vector<double>& angles, double T, double R, int seed) {

    using edge_vector = std::vector<std::pair<int, int>>;
    edge_vector result;

    std::vector<std::pair<
            edge_vector,
            uint64_t[31] /* avoid false sharing */
    > > local_edges(omp_get_max_threads());

    constexpr auto block_size = size_t{1} << 20;

    for(auto& v : local_edges)
        v.first.reserve(block_size);

    std::mutex m;
    auto flush = [&] (const edge_vector& local) {
        std::lock_guard<std::mutex> lock(m);
        result.insert(result.end(), local.cbegin(), local.cend());
    };

    auto addEdge = [&](int u, int v, int tid) {
        auto& local = local_edges[tid].first;
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
            local.reserve(block_size); // just in case
        }
    };

    auto generator = hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge);
    generator.generate(seed);

    for(const auto& v : local_edges)
        flush(v.first);

    return result;
}

} // namespace hypergirgs
