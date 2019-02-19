
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cassert>

#include <girgs/Generator.h>
#include <girgs/SpatialTree.h>
#include <girgs/ScopedTimer.h>


#include "CounterPerThread.h"

#define COUNT_EDGES_ONLY

int main(int argc, char* argv[]) {
    constexpr auto dimensions = 1;
    const auto n = 100000;
    const auto ple = 2.5;
    const auto alpha = 1.1;
    const auto deg = 10;
    const auto weightSeed = 12;
    const auto positionSeed = 130;
    const auto samplingSeed = 1400;

    const auto scaling = -36; // set to negative value to compute it

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

    auto weights = [&] {
        ScopedTimer timer("Generating weights");
        return girgs::generateWeights(n, ple, weightSeed);
    }();

    auto positions = [&] {
        ScopedTimer timer("Generating positions");
        return girgs::generatePositions(n, dimensions, positionSeed);
    }();

    if (scaling < 0.0) {
        ScopedTimer timer("Find scaling");
        auto scaler = girgs::scaleWeights(weights, deg, dimensions, alpha);
        std::cout << "Scaling: " << scaler << std::endl;
    } else {
        for(auto& w : weights)
            w *= scaling;
    }

    {
        ScopedTimer timer("Generate edges");
        girgs::makeSpatialTree<dimensions>(weights, positions, alpha, addEdge).generateEdges(samplingSeed);
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
