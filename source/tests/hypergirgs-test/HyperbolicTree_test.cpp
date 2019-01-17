
#include <algorithm>
#include <cmath>
#include <numeric>

#include <gmock/gmock.h>

#include <hypergirgs/HyperbolicTree.h>
#include <hypergirgs/Hyperbolic.h>


using namespace std;
using namespace hypergirgs;


class HyperbolicTree_test: public testing::Test
{
protected:
    const int radiiSeed = 12;
    const int angleSeed = 130;
    const int edgesSeed = 1400;
};

TEST_F(HyperbolicTree_test, testNoDuplicateEdges)
{
    const auto n = 1000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    const auto deg = 10;

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
    auto angles = hypergirgs::sampleAngles(n, angleSeed);
    auto graph = hypergirgs::generateEdges(radii, angles, T, R, edgesSeed);

    auto duplicates = 0;
    sort(graph.begin(), graph.end());
    for(int i=0; i<graph.size()-1; ++i){
        auto current = graph[i];
        auto next = graph[i+1];
        duplicates += (current == next);
    }

    ASSERT_EQ(duplicates, 0);
}


TEST_F(HyperbolicTree_test, testThresholdModel)
{
    const auto n = 1000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0;
    const auto deg = 10;

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
    auto angles = hypergirgs::sampleAngles(n, angleSeed);
    auto graph = hypergirgs::generateEdges(radii, angles, T, R, edgesSeed);

    // convert to adjacency list
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
            auto should_be_present = hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]) < R;
            ASSERT_EQ(edge_present, should_be_present);
        }
    }
}

TEST_F(HyperbolicTree_test, testGeneralModel)
{
    const auto n = 1000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto T = 0.5;
    const auto deg = 10;

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
    auto angles = hypergirgs::sampleAngles(n, angleSeed);
    auto graph = hypergirgs::generateEdges(radii, angles, T, R, edgesSeed);

    auto num_edges = graph.size();
    auto num_desired = 0.5*deg*n;
    auto rigor = 0.9;
    ASSERT_LE(rigor * num_edges, num_desired);
    ASSERT_LE(rigor * num_desired, num_edges);
}


TEST_F(HyperbolicTree_test, testReproducible)
{
    const auto n = 1000;
    const auto alpha = 0.75; // ple = 2*alpha+1
    const auto Ts = {0.0, 0.5};
    const auto deg = 10;

    for(auto T : Ts) {
        auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
        auto radii = hypergirgs::sampleRadii(n, alpha, R, radiiSeed);
        auto angles = hypergirgs::sampleAngles(n, angleSeed);
        auto edges1 = hypergirgs::generateEdges(radii, angles, T, R, edgesSeed);
        auto edges2 = hypergirgs::generateEdges(radii, angles, T, R, edgesSeed);
        sort(edges1.begin(), edges1.end());
        sort(edges2.begin(), edges2.end());
        ASSERT_EQ(edges1, edges2);
    }
}
