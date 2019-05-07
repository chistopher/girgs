
#include <hypergirgs/Generator.h>

#include <random>
#include <fstream>
#include <cmath>
#include <mutex>

#include <omp.h>

#include <hypergirgs/HyperbolicTree.h>


namespace hypergirgs {


double calculateRadius(int n, double alpha, double T, double deg) {
    return 2 * log(n * 2 * alpha * alpha * (T == 0 ? 1 / PI : T / sin(PI * T)) /
                   (deg * (alpha - 0.5) * (alpha - 0.5)));
}

///////////////////////////////// COPIED FROM NETWORKIT ///////////////////////
static double getExpectedDegree(double n, double alpha, double R) {
    double gamma = 2*alpha+1;
    double xi = (gamma-1)/(gamma-2);
    double firstSumTerm = exp(-R/2);
    double secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((PI/4)*pow((1/alpha),2)-(PI-1)*(1/alpha)+(PI-2))-1);
    double expectedDegree = (2/PI)*xi*xi*n*(firstSumTerm + secondSumTerm);
    return expectedDegree;
}

static double searchTargetRadiusForColdGraphs(double n, double k, double alpha, double epsilon) {
    double gamma = 2*alpha+1;
    double xiInv = ((gamma-2)/(gamma-1));
    double v = k * (PI/2)*xiInv*xiInv;
    double currentR = 2*log(n / v);
    double lowerBound = currentR/2;
    double upperBound = currentR*2;
    assert(getExpectedDegree(n, alpha, lowerBound) > k);
    assert(getExpectedDegree(n, alpha, upperBound) < k);
    do {
        currentR = (lowerBound + upperBound)/2;
        double currentK = getExpectedDegree(n, alpha, currentR);
        // std::cout << "n: " << n << " k: " << k << " alpha: " << alpha << " curK: " << currentK << " R: " << currentR << std::endl;
        if (currentK < k) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while (std::abs(getExpectedDegree(n, alpha, currentR) - k) > epsilon );
    return currentR;
}

static double getTargetRadius(double n, double m, double alpha=1, double T=0, double epsilon = 0.01) {
    double result;
    double plexp = 2*alpha+1;
    double targetAvgDegree = (m/n)*2;
    double xiInv = ((plexp-2)/(plexp-1));
    if (T == 0) {
        double v = targetAvgDegree * (PI/2)*xiInv*xiInv;
        result = 2*log(n / v);
        result = searchTargetRadiusForColdGraphs(n, targetAvgDegree, alpha, epsilon);
    } else {
        double beta = 1/T;
        if (T < 1){//cold regime
            double Iinv = ((beta/PI)*sin(PI/beta));
            double v = (targetAvgDegree*Iinv)*(PI/2)*xiInv*xiInv;
            result = 2*log(n / v);
        } else {//hot regime
            double v = targetAvgDegree*(1-beta)*pow((PI/2), beta)*xiInv*xiInv;
            result = 2*log(n/v)/beta;
        }
    }
    return result;
}
////////////////////////////////// END NETWORKIT COPY ////////////////////

double calculateRadiusLikeNetworKit(int n, double alpha, double T, double deg) {
    return getTargetRadius(n, 0.5*deg*n, alpha, T);
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
        auto adist = std::uniform_real_distribution<>(0, 2*PI);
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
