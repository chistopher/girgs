
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
    const auto alpha = numeric_limits<double>::infinity();
    const auto ple = -2.8;

    Generator generator;
    generator.setWeights(n, ple, seed);
    auto weights = generator.weights();
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){

        generator.setPositions(n, d, seed+d);
        generator.generateTreshold();

        // check that there is an edge if and only if the condition in the paper holds: dist < c*(w1w2/W)^-d
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){
                auto& a = generator.graph()[j];
                auto& b = generator.graph()[i];

                auto dist = distance(a.coord, b.coord);
                auto w = std::pow(a.weight * b.weight / W, 1.0/d);

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
    const auto alpha = 2.5;
    const auto ple = -2.5;

    auto generator = Generator();
    generator.setWeights(n, ple, seed);
    auto weights = generator.weights();
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){
        // check that the number of generated edges is close to the expected value

        // 1) generator
        generator.setPositions(n, d, seed+d);
        generator.generate(alpha, seed+d);

        // 2) quadratic sanity check
        auto expectedEdges = vector<double>(n, 0.0);
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){
                auto& a = generator.graph()[j];
                auto& b = generator.graph()[i];

                auto dist = std::pow(distance(a.coord, b.coord), d);
                auto w = a.weight * b.weight / W;

                auto prob = std::min(std::pow(w/dist, alpha), 1.0);
                expectedEdges[i] += prob;
                expectedEdges[j] += prob;
            }
        }

        auto total_expected = accumulate(expectedEdges.begin(), expectedEdges.end(), 0.0);
        auto total_actual = accumulate(generator.graph().begin(), generator.graph().end(), 0.0, [](double sum, const Node& node){
            return sum + node.edges.size();
        });

        auto rigor = 0.99;
        test(rigor * total_expected < total_actual);
        test(rigor * total_actual < total_expected);
    }
}


void testCompleteGraph(int seed) {

    const auto n = 100;
    const auto alpha = 0.0; // each edge prob will be 100% now
    const auto ple = -2.5;

    auto generator = Generator();
    generator.setWeights(n, ple, seed);

    for(auto d=1u; d<5; ++d) {

        generator.setPositions(n, d, seed+d);
        generator.generate(alpha, seed+d);

        // check that each node is connected to all other nodes
        for(auto& node : generator.graph()) {
            test(node.edges.size() == n-1);
            for(auto& other : generator.graph())
                if(node.index != other.index)
                    test(find(node.edges.begin(), node.edges.end(), &other) != node.edges.end());
        }
    }
}



// samples all edges by threshold model: dist(i,j) < c*(wiwj/W)^(1/d)
double edgesInQuadraticSampling(const std::vector<double>& w, const vector<vector<double>>& pos, double c) {
    auto n = w.size();
    auto d = pos.front().size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto edges = 0.0;
    for(int i=0; i<n; ++i)
        for(int j=i+1; j<n; ++j)
            if(distance(pos[i], pos[j]) < c*std::pow(w[i] * w[j] / W, 1.0/d))
                edges += 2; // both endpoints get an edge
    return edges;
}

void testThresholdEstimation(int seed) {

    auto n = 100;
    auto PLE = -2.5;
    auto alpha = numeric_limits<double>::infinity();
    auto weightSeed = seed;
    auto positionSeed = seed;

    auto desired_avg = 10;
    auto runs = 50;

    Generator generator;
    generator.setWeights(n, PLE, weightSeed);
    auto weights = generator.weights();

    // do the tests for all dimensions < 5
    for(auto d = 1; d<5; ++d) {

        // estimate scaling for current dimension
        generator.setWeights(weights); // reset weights
        auto scaling = generator.scaleWeights(desired_avg, d, alpha);
        auto estimated_c = pow(scaling, 1.0/d);

        // observed avg with estimated c (over multiple runs with different positions)
        auto observed_avg = 0.0;
        for(int i = 0; i<runs; ++i) {

            // try GIRGS generator and quadratic sampling
            generator.setPositions(n, d, positionSeed+i);
            generator.generateTreshold();

            auto avg1 = generator.avg_degree();
            auto avg2 = edgesInQuadraticSampling(weights, generator.positions(), estimated_c) / n;

            // generator must yield same results as quadratic sampling
            test(avg1 == avg2);
            observed_avg += avg1;
        }
        observed_avg /= runs;

        // test the goodness of the estimation for weight scaling
        test(abs(desired_avg - observed_avg) < 0.1);
    }
}


void testWeightSampling(int seed) {

    auto n = 10000;
    auto ple = -2.1;
    int runs = 10;

    Generator g;

    for(int i=0; i<runs; ++i){
        g.setWeights(n, ple, seed+i);
        auto weights = g.weights();

        for(auto each : weights) {
            test(each >= 1.0);
            test(each < n);
        }
        auto max_weight = *max_element(weights.begin(), weights.end());
        test(max_weight > n/10);
    }
}



int main(int argc, char* argv[]) {

    const auto seed = 1337;

    testThresholdModel(seed);
    testGeneralModel(seed);
    testCompleteGraph(seed);
    testThresholdEstimation(seed);
    testWeightSampling(seed);

    cout << "all tests passed." << endl;
    return 0;
}
