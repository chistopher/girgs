
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cassert>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>


using namespace std;
using namespace chrono;

int main(int argc, char* argv[]) {
    const auto n = 1000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    const auto deg = 10;
    const auto radiiSeed = 12;
    const auto angleSeed = 130;
    const auto edgesSeed = 1400;

    std::vector<std::pair<int,int>> graph;
    auto addEdge = [&graph] (int u, int v, int tid) {
        assert(tid == 0);
        graph.emplace_back(u,v);
    };

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
    auto angles = hypergirgs::sampleAngles(n, angleSeed);
    auto generator = hypergirgs::makeHyperbolicTree(radii, angles, T, R, addEdge);

    // measure
    auto start = std::chrono::high_resolution_clock::now();
    generator.generate(edgesSeed);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}
