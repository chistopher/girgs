
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cassert>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/ScopedTimer.h>

#include "CounterPerThread.h"

#define COUNT_EDGES_ONLY

int main(int argc, char* argv[]) {
    const auto n = 1000000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    const auto deg = 10;
    const auto radiiSeed = 12;
    const auto angleSeed = 130;
    const auto edgesSeed = 1400;

#ifdef COUNT_EDGES_ONLY
    // only count edges
    CounterPerThread<uint64_t> num_edges;

    auto addEdge = [&num_edges] (int /* u */, int /* v */, int tid) {
        num_edges.add(tid);
    };

    std::cout << "!! Count edges only !!\n";
#else
    // actually store graph
    std::vector<std::pair<int,int>> graph;
    graph.reserve(n*deg/2);
    auto addEdge = [&graph] (int u, int v, int tid) {
        assert(tid == 0);
        graph.emplace_back(u,v);
    };
#endif

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);

    double time_points, time_preprocess, time_sample;

    // Sample points
    std::vector<double> radii, angles;
    {
        ScopedTimer timer("Generate points", time_points);
        radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
        angles = hypergirgs::sampleAngles(n, angleSeed);
    }

    // Preprocess
    auto generator = [&] {
        ScopedTimer timer("Preprocess", time_preprocess);
        return hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge, true);
    }();

    // Generate edges
    {
        ScopedTimer timer("Generate edges", time_sample);
        generator.generate(edgesSeed);
    }

    const auto total_edges =
#ifdef COUNT_EDGES_ONLY
        num_edges.total();
#else
        graph.size();
#endif

    std::cout << "Number of edges generated: " << total_edges << "\n"
                 "Average degree: " << (2.0 * total_edges / n) << "\n";
}
