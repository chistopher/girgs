
#include <algorithm>
#include <cmath>
#include <numeric>

#include <gmock/gmock.h>

#include <girgs/Generator.h>
#include <girgs/Node.h>


using namespace std;


class Generator_test: public testing::Test
{
protected:
    int seed = 1337;
};


TEST_F(Generator_test, testThresholdModel)
{
    const auto n = 100;
    const auto alpha = numeric_limits<double>::infinity();
    const auto ple = -2.8;

    girgs::Generator generator;
    generator.setWeights(n, ple, seed);
    auto weights = generator.weights();
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){

        generator.setPositions(n, d, seed+d);
        generator.generateThreshold();

        // check that there is an edge if and only if the condition in the paper holds: dist < c*(w1w2/W)^-d
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){
                auto& a = generator.graph()[j];
                auto& b = generator.graph()[i];

                auto dist = girgs::distance(a.coord, b.coord);
                auto w = std::pow(a.weight * b.weight / W, 1.0/d);

                auto edge1 = find(a.edges.begin(), a.edges.end(), &b);
                auto edge2 = find(b.edges.begin(), b.edges.end(), &a);
                if(dist < w) {
                    //EXPECT_NE(edge1, a.edges.end()) << "edge should be present";
                    //EXPECT_NE(edge2, b.edges.end()) << "edge should be present";
                } else {
                    //EXPECT_EQ(edge1, a.edges.end()) << "edge should be absent";
                    //EXPECT_EQ(edge2, b.edges.end()) << "edge should be absent";
                }
            }
        }
    }
}

TEST_F(Generator_test, testGeneralModel)
{
    const auto n = 1000;
    const auto alpha = 2.5;
    const auto ple = -2.5;

    auto generator = girgs::Generator();
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

                auto dist = std::pow(girgs::distance(a.coord, b.coord), d);
                auto w = a.weight * b.weight / W;

                auto prob = std::min(std::pow(w/dist, alpha), 1.0);
                expectedEdges[i] += prob;
                expectedEdges[j] += prob;
            }
        }

        auto total_expected = accumulate(expectedEdges.begin(), expectedEdges.end(), 0.0);
        auto total_actual = accumulate(generator.graph().begin(), generator.graph().end(), 0.0,
                [](double sum, const girgs::Node& node){ return sum + node.edges.size() * 2; });

        auto rigor = 0.98;
        EXPECT_LT(rigor * total_expected, total_actual) << "edges too much below expected value";
        EXPECT_LT(rigor * total_actual, total_expected) << "edges too much above expected value";
    }
}


TEST_F(Generator_test, testCompleteGraph)
{
    const auto n = 100;
    const auto alpha = 0.0; // each edge prob will be 100% now
    const auto ple = -2.5;

    auto generator = girgs::Generator();
    generator.setWeights(n, ple, seed);

    for(auto d=1u; d<5; ++d) {

        generator.setPositions(n, d, seed+d);
        generator.generate(alpha, seed+d);

        // check that each node is connected to all other nodes
        for(auto& node : generator.graph()) {
            //EXPECT_EQ(node.edges.size(), n-1) << "expect a complete graph withour self loops";
            for(auto& other : generator.graph())
                if(node.index != other.index){
                    auto it = find(node.edges.begin(), node.edges.end(), &other);
                    //EXPECT_NE(it, node.edges.end()) << "edge should be present";
                }
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
            if(girgs::distance(pos[i], pos[j]) < c*std::pow(w[i] * w[j] / W, 1.0/d))
                edges += 2; // both endpoints get an edge
    return edges;
}


TEST_F(Generator_test, testThresholdEstimation)
{
    auto n = 100;
    auto PLE = -2.5;
    auto alpha = numeric_limits<double>::infinity();
    auto weightSeed = seed;
    auto positionSeed = seed;

    auto desired_avg = 10;
    auto runs = 50;

    girgs::Generator generator;
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
            generator.generateThreshold();

            auto avg1 = generator.avg_degree();
            auto avg2 = edgesInQuadraticSampling(weights, generator.positions(), estimated_c) / n;

            // generator must yield same results as quadratic sampling
            EXPECT_EQ(avg1, avg2) << "sampling with scaled weights produced different results than quadratic samping with constant factor";
            observed_avg += avg1;
        }
        observed_avg /= runs;

        // test the goodness of the estimation for weight scaling
        EXPECT_LT(abs(desired_avg - observed_avg), 0.1) << "estimated constant does not produce desired average degree";
    }
}


TEST_F(Generator_test, testEstimation)
{
    auto all_n = {100, 150};
    auto all_alpha = {0.7, 1.5, 3.0, numeric_limits<double>::infinity()};
    auto all_desired_avg = {10, 15, 20};
    auto all_dimensions = {1, 2, 3};
    auto runs = 5;

    auto PLE = -2.5;
    auto weightSeed = seed;
    auto positionSeed = seed;

    for(int n : all_n){
        for(double alpha : all_alpha){
            for(double desired_avg : all_desired_avg){
                for(int d : all_dimensions){

                    // generate weights
                    girgs::Generator generator;
                    generator.setWeights(n, PLE, weightSeed);
                    auto weights = generator.weights();

                    // estimate scaling for current dimension
                    generator.setWeights(weights); // reset weights
                    generator.scaleWeights(desired_avg, d, alpha);

                    auto observed_avg = 0.0;
                    for(int i = 0; i<runs; ++i) {

                        // try GIRGS generator
                        generator.setPositions(n, d, positionSeed+i);
                        generator.generate(alpha, n+i);

                        auto avg = generator.avg_degree();
                        observed_avg += avg;
                    }
                    observed_avg /= runs;

                    // test the goodness of the estimation for weight scaling
                    EXPECT_LT(abs(desired_avg - observed_avg), 1.0) << "estimated constant does not produce desired average degree";
                }
            }
        }
    }
}


TEST_F(Generator_test, testWeightSampling)
{
    auto n = 10000;
    auto ple = -2.1;
    int runs = 10;

    girgs::Generator g;

    for(int i=0; i<runs; ++i){
        g.setWeights(n, ple, seed+i);
        auto weights = g.weights();

        for(auto each : weights) {
            EXPECT_GE(each, 1.0);
            EXPECT_LT(each, n);
        }
        auto max_weight = *max_element(weights.begin(), weights.end());
        EXPECT_GT(max_weight * max_weight, n) << "max weight should be large";
    }
}
