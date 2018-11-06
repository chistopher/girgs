
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

#include <Generator.h>


using namespace std;

void test(bool cond){
    if(!cond){
        cout << "assertion failed" << endl;
        exit(1);
    }
}


void testThresholdModel(int seed) {

    const auto n = 100;
    const auto c = 1.0;
    const auto alpha = numeric_limits<double>::infinity();

    auto weights = generateWeights(n, -2.5, seed);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);
    auto generator = Generator();

    for(auto d=1u; d<5; ++d){

        generator.generateGIRG(d, weights, alpha, c, seed);

        // check that there is an edge if and only if the condition in the paper holds: dist < c*(w1w2/W)^-d
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){
                auto& a = generator.graph()[j];
                auto& b = generator.graph()[i];

                auto dist = distance(a.coord, b.coord);
                auto w = c*std::pow(a.weight * b.weight / W, 1.0/d);

                auto edge1 = find(a.edges.begin(), a.edges.end(), &b);
                auto edge2 = find(b.edges.begin(), b.edges.end(), &a);
                if(dist < w) {
                    test(edge1 != a.edges.end());
                    test(edge2 != b.edges.end());
                } else {
                    test(edge1 == a.edges.end());
                    test(edge2 == b.edges.end());
                }
            }
        }
    }
}


void testGeneralModel(int seed) {

    const auto n = 1000;
    const auto c = 0.75;
    const auto alpha = 1.5;

    auto weights = generateWeights(n, -2.5, seed);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);
    auto generator = Generator();

    for(auto d=1u; d<5; ++d){

        generator.generateGIRG(d, weights, alpha, c, seed);

        // check that the number of generated edges is close to the expected value
        auto expectedEdges = vector<double>(n, 0.0);
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){
                auto& a = generator.graph()[j];
                auto& b = generator.graph()[i];

                auto dist = std::pow(distance(a.coord, b.coord), d);
                auto w = a.weight * b.weight / W;

                auto prob = std::min(c*std::pow(w/dist, alpha), 1.0);
                expectedEdges[i] += prob;
                expectedEdges[j] += prob;
            }
        }

        auto total_expected = accumulate(expectedEdges.begin(), expectedEdges.end(), 0.0);
        auto total_actual = 0.0;
        for(auto& node : generator.graph()) {
            auto expected = expectedEdges[node.index];
            auto actual = node.edges.size();
            total_actual += actual;
        }

        auto rigor = 0.99;
        test(rigor * total_expected < total_actual);
        test(rigor * total_actual < total_expected);
    }
}


void testCompleteGraph(int seed) {

    const auto n = 100;
    const auto c = 1.0;
    const auto alpha = 0.0; // each edge prob will be 100% now

    auto weights = generateWeights(n, -2.5, seed);
    auto generator = Generator();

    for(auto d=1u; d<5; ++d) {

        generator.generateGIRG(d, weights, alpha, c, seed);

        // check that each node is connected to all other nodes
        for(auto& node : generator.graph()) {
            test(node.edges.size() == n-1);

            for(auto& other : generator.graph())
                if(node.index != other.index)
                    test(find(node.edges.begin(), node.edges.end(), &other) != node.edges.end());
        }
    }
}




int main(int argc, char* argv[]) {

    const auto seed = 1337;

    testThresholdModel(seed);
    testGeneralModel(seed);
    testCompleteGraph(seed);

    cout << "all tests passed." << endl;
    return 0;
}
