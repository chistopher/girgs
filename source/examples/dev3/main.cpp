
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

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
    auto angles = hypergirgs::sampleAngles(n, angleSeed);
    auto generator = hypergirgs::HyperbolicTree(radii, angles, T, R);

    // measure
    auto start = std::chrono::high_resolution_clock::now();
    auto graph = generator.generate(edgesSeed);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    // test no duplicate edges
    auto duplicates = 0;
    sort(graph.begin(), graph.end());
    for(int i=0; i<graph.size()-1; ++i){
        auto current = graph[i];
        auto next = graph[i+1];
        duplicates += (current == next);
    }

    cout << "edges:      " << graph.size() << endl;
    cout << "avg deg:    " << 2.0* graph.size() / n << endl;
    cout << "duplicates: " << duplicates << endl;

    // preliminary tests
    cout << "computing stats ...\t\t" << flush;
    auto t6 = high_resolution_clock::now();

    // init accumulators
    long long missing_edges = 0;
    long long false_edges = 0;
    auto missing_deviation = 0.0;
    auto false_deviation = 0.0;

    // normalize adj matrix that edges point from smaller index to larger
    auto adj_list = vector<vector<int>>(n);
    for(auto edge : graph) {
            auto mm = minmax(edge.first, edge.second);
            adj_list[mm.first].push_back(mm.second);
    }

    // check all point pairs
    auto current_neighs = std::vector<char>(n);
    for(int i=0; i<n; ++i){
        // get one row in adj matrix
        std::fill(current_neighs.begin(), current_neighs.end(), false);
        for(auto neigh : adj_list[i])
            current_neighs[neigh] = true;
        // check that all edges of i are correct in hyperbolic space
        for(int j=i+1; j<n; ++j){
            auto edge_present = current_neighs[j];
            auto dist = hypergirgs::hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]);
            if(dist < R) {
                if(!edge_present) {
                    missing_deviation += R - dist;
                    missing_edges++;
                }
            } else {
                if(edge_present) {
                    false_deviation += dist - R;
                    false_edges++;
                }
            }
        }
    }

    auto t7 = high_resolution_clock::now();
    auto total_edges = graph.size();
    cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    cout << '\n';
    cout << "total edges       " << total_edges << '\n';
    cout << '\n';
    cout << "missing edges     " << missing_edges << '\n';
    cout << "missing deviation " << missing_deviation << '\n';
    cout << "average mistake   " << missing_deviation / missing_edges << '\n';
    cout << '\n';
    cout << "false edges       " << false_edges << '\n';
    cout << "false deviation   " << false_deviation << '\n';
    cout << "average mistake   " << false_deviation / false_edges << '\n';
    cout << '\n';
    cout << "(false+missing)/total " << static_cast<double>(false_edges+missing_edges)/total_edges << '\n';
}
