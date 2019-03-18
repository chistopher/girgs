#include <iostream>
#include <sstream>
#include <iomanip>

#include <algorithm>
#include <cassert>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/ScopedTimer.h>

#include "CounterPerThread.h"

#ifdef __unix__
#include "unistd.h"

static std::string hostname() {
    char tmp[128];
    if (gethostname(tmp, 128))
        return "n/a";
    return {tmp};
}

#else
static std::string hostname() {
    return "n/a";
}
#endif

/////////////////////////////////////////////////////////////
// copied from NetworKit
static double getExpectedDegree(double n, double alpha, double R) {
    double gamma = 2*alpha+1;
    double xi = (gamma-1)/(gamma-2);
    double firstSumTerm = exp(-R/2);
    double secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((M_PI/4)*pow((1/alpha),2)-(M_PI-1)*(1/alpha)+(M_PI-2))-1);
    double expectedDegree = (2/M_PI)*xi*xi*n*(firstSumTerm + secondSumTerm);
    return expectedDegree;
}

static double searchTargetRadiusForColdGraphs(double n, double k, double alpha, double epsilon) {
    double gamma = 2*alpha+1;
    double xiInv = ((gamma-2)/(gamma-1));
    double v = k * (M_PI/2)*xiInv*xiInv;
    double currentR = 2*log(n / v);
    double lowerBound = currentR/2;
    double upperBound = currentR*2;
    assert(getExpectedDegree(n, alpha, lowerBound) > k);
    assert(getExpectedDegree(n, alpha, upperBound) < k);
    do {
        currentR = (lowerBound + upperBound)/2;
        double currentK = getExpectedDegree(n, alpha, currentR);
        std::cout << "n: " << n << " k: " << k << " alpha: " << alpha << " curK: " << currentK << " R: " << currentR << std::endl;
        if (currentK < k) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while (abs(getExpectedDegree(n, alpha, currentR) - k) > epsilon );
    return currentR;
}

static double getTargetRadius(double n, double m, double alpha=1, double T=0, double epsilon = 0.01) {
    double result;
    double plexp = 2*alpha+1;
    double targetAvgDegree = (m/n)*2;
    double xiInv = ((plexp-2)/(plexp-1));
    if (T == 0) {
        double v = targetAvgDegree * (M_PI/2)*xiInv*xiInv;
        result = 2*log(n / v);
        result = searchTargetRadiusForColdGraphs(n, targetAvgDegree, alpha, epsilon);
    } else {
        double beta = 1/T;
        if (T < 1){//cold regime
            double Iinv = ((beta/M_PI)*sin(M_PI/beta));
            double v = (targetAvgDegree*Iinv)*(M_PI/2)*xiInv*xiInv;
            result = 2*log(n / v);
        } else {//hot regime
            double v = targetAvgDegree*(1-beta)*pow((M_PI/2), beta)*xiInv*xiInv;
            result = 2*log(n/v)/beta;
        }
    }
    return result;
}
/////////////////////////////////////////////////////////////


double benchmark(std::ostream& os, const std::string& host, unsigned int iter, unsigned int n, unsigned int avgDeg, double alpha, double T, unsigned int seed = 0) {
    CounterPerThread<uint64_t> counter_num_edges;

    double time_total, time_points, time_preprocess, time_sample;

    const auto R = getTargetRadius(n, avgDeg * n / 2.0, alpha, T);
    {
        ScopedTimer tot_timer("Total", time_total);
        //auto R = hypergirgs::calculateRadius(n, alpha, T, avgDeg);

        // Sample points
        std::vector<double> radii, angles;
        {
            ScopedTimer timer("Generate points", time_points);
            std::tie(radii, angles) = hypergirgs::sampleRadiiAndAngles(n, alpha, R, seed, true);
        }

        // Preprocess
        auto addEdge = [&counter_num_edges] (int, int, int tid) {counter_num_edges.add(tid);};
        hypergirgs::HyperbolicTree<decltype(addEdge)> generator = [&] {
            ScopedTimer timer("Preprocess", time_preprocess);
            return hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge, true);
        }();

        // Generate edges
        {
            ScopedTimer timer("Generate edges", time_sample);
            generator.generate(seed+12345678);
        }
    }

    const auto num_edges = counter_num_edges.total();

    // Logging
    {
        std::stringstream ss;
        ss << "[CSV]"
           << "hypergirgs,"
           << host << ","
           << iter << ","
           << seed << ","
           << n << ","
           << avgDeg << ","
           << alpha << ","
           << (2*alpha + 1) << ","
           << T << ","
           << R << ","
           << time_total      << ","
           << time_points     << ","
           << time_preprocess << ","
           << time_sample << ","
           << num_edges << ","
           << (2.0 * num_edges / n);

        os << ss.str() << std::endl;
    }

    return time_total;
}

int main(int argc, char* argv[]) {
    unsigned int seed = 0;

    const unsigned n0 = 1e4;
    const unsigned nMax = 1e8;
    const unsigned steps_per_dec = 3;
    const double timeout = 100 * 1e3; // ms


    // Print Header
    std::cerr <<
        "[CSV]"
          "algo,"
          "host,"
          "iter,"
          "seed,"
          "n,"
          "avgDeg,"
          "alpha,"
          "PLE,"
          "T,"
          "R,"
          "TimeTotal,"
          "TimePoints,"
          "TimePrepro,"
          "TimeEdges,"
          "GenNumEdge,"
          "GenAvgDeg\n";

    const auto host = hostname();

    for(int iter = 0; iter < 5; iter++) {
        for (const double T : {0.0, 0.5, 0.9}) {
            for (const double ple : {2.2, 3.0}) {
                unsigned int skip_n = nMax + 1;
                const auto alpha = (ple - 1.0) / 2.0;

                for (const int avgDeg : {10, 100, 1000}) {
                    int ni = 1;

                    for (auto n = n0; n <= nMax; n = n0 * std::pow(10.0, 1.0 * ni / steps_per_dec), ++ni) {
                        if (avgDeg * 20 > n)
                            continue;

                        std::cout << "\033[1miter=" << iter << ", T=" << T << ", PLE=" << ple << ", n=" << n << ", avgDeg=" << avgDeg
                                  << "\033[21m\n";

                        double time;
                        if (n < skip_n) {
                            // if last (smaller) problem took too long, we skip this one
                            // we do not use break to make sure seed stays consistent
                            time = benchmark(std::cerr, host, iter, n, avgDeg, alpha, T, seed);
                            if (time > timeout) {
                                skip_n = n;
                                std::cout << " took too long\n";
                            }
                        } else {
                            std::cout << " skip_n = " << skip_n << "\n";
                        }

                        seed += 10;
                    }
                }
            }
        }
    }

    return 0;
}
