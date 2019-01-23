#include <iostream>
#include <sstream>
#include <iomanip>

#include <algorithm>
#include <cassert>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>

#include <CounterPerThread.h>
#include <ScopedTimer.h>

void benchmark(std::ostream& os, unsigned int n, unsigned int avgDeg, double alpha, double T, unsigned int seed = 0) {
    CounterPerThread<uint64_t> counter_num_edges;

    double time_total, time_points, time_preprocess, time_sample;
    {
        ScopedTimer tot_timer("Total", time_total);
        auto R = hypergirgs::calculateRadius(n, alpha, T, avgDeg);

        // Sample points
        std::vector<double> radii, angles;
        {
            ScopedTimer timer("Generate points", time_points);
            radii = hypergirgs::sampleRadii(n, alpha, R, seed++);
            angles = hypergirgs::sampleAngles(n, seed++);
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
            generator.generate(seed++);
        }
    }

    const auto num_edges = counter_num_edges.total();

    // Logging
    {
        std::stringstream ss;
        ss << "[DATA] "
           << std::setw(11) << "SimplePP,"
           << std::setw(10) << n << ","
           << std::setw(10) << avgDeg << ","
           << std::setw(10) << alpha << ","
           << std::setw(10) << T << ","
           << std::setw(10) << time_total      << ","
           << std::setw(10) << time_points     << ","
           << std::setw(10) << time_preprocess << ","
           << std::setw(10) << time_sample << ","
           << std::setw(10) << num_edges << ","
           << std::setw(10) << (2.0 * num_edges / n);

        os << ss.str() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    unsigned int seed = 0;

    // Print Header
    {
        std::stringstream ss;
        ss << "[DATA] "
           << std::setw(11) << "algo,"
           << std::setw(11) << "n,"
           << std::setw(11) << "avgDeg,"
           << std::setw(11) << "alpha,"
           << std::setw(11) << "T,"
           << std::setw(11) << "TimeTotal,"
           << std::setw(11) << "TimePoints,"
           << std::setw(11) << "TimePrepro,"
           << std::setw(11) << "TimeEdges,"
           << std::setw(11) << "GenNumEdge,"
           << std::setw(10) << "GenAvgDeg";


        std::cout << ss.str() << std::endl;
    }


    for(unsigned int iter = 0; iter != 5; ++iter) {
        for (unsigned int n = 1024; n < (2 << 21); n *= 2) {
            for (unsigned avgDeg : {10, 100, 1000}) {
                if (avgDeg * 10 > n) continue;

                std::clog << "iter=" << iter << ", n=" << n << ", avgDeg=" << avgDeg << "\n";

                benchmark(std::cout, n, avgDeg, alpha, T, seed);
                seed += 10;
            }
        }
    }

    return 0;
}
