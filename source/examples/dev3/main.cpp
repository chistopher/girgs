
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cassert>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>


using namespace std;
using namespace chrono;

int main(int argc, char* argv[]) {
    const auto n = 1000000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    const auto deg = 10;
    const auto radiiSeed = 12;
    const auto angleSeed = 130;
    const auto edgesSeed = 1400;

    using counter_with_padding = std::pair<uint64_t, char[64]>;
    std::vector<counter_with_padding> num_edges_per_thread(omp_get_max_threads());
    for(auto& x: num_edges_per_thread) x.first = 0;

    auto addEdge = [&num_edges_per_thread] (int u, int v, int tid) {
        num_edges_per_thread[tid].first++;
    std::vector<std::pair<int,int>> graph;
    graph.reserve(n*deg/2);
    auto addEdge = [&graph] (int u, int v, int tid) {
        assert(tid == 0);
        graph.emplace_back(u,v);
    };

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
        return hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge);
    }();

    // Generate edges
    {
        ScopedTimer timer("Generate edges", time_sample);
        generator.generate(edgesSeed);
    }


}
