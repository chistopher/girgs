
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

template <bool Radii, bool Angles>
static std::pair<std::vector<double>, std::vector<double>> sampleRadiiAndAnglesHelper(
    const int n, const double alpha, const double R, const int seed, const bool parallel
) {
    static_assert(Radii || Angles, "At least one output is required");

    std::vector<double> radii(n * Radii);
    std::vector<double> angles(n * Angles);

    constexpr auto kMinChunkSize = 10000;
    const auto threads = parallel ? std::min<int>(omp_get_max_threads(), (n + kMinChunkSize - 1) / kMinChunkSize) : 1;

    const auto invalpha = 1.0 / alpha;
    #pragma omp parallel num_threads(threads)
    {
        auto gen = hypergirgs::default_random_engine(seed >= 0 ? seed + omp_get_thread_num() : std::random_device()());
        auto adist = std::uniform_real_distribution<>(0, 2*M_PI);
        auto rdist = std::uniform_real_distribution<>(std::nextafter(1.0, 2.0), std::cosh(alpha * R));

        // warm-up generator
        constexpr int kMinWarmup = 1000;
        for (int i = 0; i < std::max<int>(n / threads / 5, kMinWarmup); ++i)
            gen();

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            if (Angles)
                angles[i] = adist(gen);

            if (Radii)
                radii[i] = acosh(rdist(gen)) * invalpha;
        }
    }

    return {radii, angles};
}


std::vector<double> sampleRadii(int n, double alpha, double R, int seed, bool parallel) {
    return sampleRadiiAndAnglesHelper<true, false>(n, alpha, R, seed, parallel).first;
}

std::vector<double> sampleAngles(int n, int seed, bool parallel) {
    return sampleRadiiAndAnglesHelper<false, true>(n, /*unused*/1.0, /*unused*/10.0, seed, parallel).second;

}

std::pair<std::vector<double>, std::vector<double>> sampleRadiiAndAngles(int n, double alpha, double R, int seed, bool parallel) {
    return sampleRadiiAndAnglesHelper<true, true>(n, alpha, R, seed, parallel);
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
