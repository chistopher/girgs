
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <Generator.h>


using namespace std;

void test(bool cond){
    if(!cond){
        cout << "assertion failed" << endl;
        exit(1);
    }
}


double estimate_scaling(const std::vector<double>& weights, int desired_avg_degree, int dimension) {

    // currently only for threshold model!!!!

    // compute some constant stuff
    auto max_weight = *max_element(weights.begin(), weights.end());
    auto W = 0.0, sq_W = 0.0;
    for(auto each : weights){
        W += each;
        sq_W += each*each;
    }

    // my function to do the exponential search on
    auto f = [W, sq_W, &weights, dimension, n(weights.size()), max_weight](double c) {
        // compute rich club
        vector<double> rich_club;
        for(auto weight : weights)
            if(std::pow(2*c,dimension) * (weight*max_weight/W) > 1.0)
                rich_club.push_back(weight);
        sort(rich_club.begin(), rich_club.end(), greater<>());
        // compute overestimation
        auto overestimation = pow(2, dimension) * pow(c, dimension) / n * (W - sq_W/W);
        // subtract error
        auto error = 0.0;
        for(int i = 0; i<rich_club.size(); ++i)
            for(int j = 0; j<rich_club.size(); ++j) {
                if(i==j) continue;
                auto w1 = rich_club[i];
                auto w2 = rich_club[j];
                auto e = max( std::pow(2*c,dimension)*(w1*w2/W)-1.0, 0.0);
                error += e;
                if(e <= 0) break;
            }
        return overestimation - error/n;
    };

    // do exponential search on avg_degree function
    auto upper = pow( weights.size()*desired_avg_degree / pow(2.0,dimension) / (W-sq_W/W), 1.0/dimension);
    auto lower = upper / 2.0;

    // scale interval up if necessary
    while(f(upper) < desired_avg_degree){
        lower = upper;
        upper *= 2;
    }

    // scale interval down if necessary
    while(f(lower) > desired_avg_degree){
        upper = lower;
        lower /= 2;
    }

    // do binary search
    auto mid = f((upper+lower)/2);
    while(abs(mid - desired_avg_degree) > 0.02) {
        if(mid < desired_avg_degree)
            lower = (upper+lower)/2;
        else
            upper = (upper+lower)/2;
        mid = f((upper+lower)/2);
    }

    return (upper+lower)/2;
}


double quadraticSampling(const std::vector<double>& w, double c, int d, const vector<vector<double>>& pos) {

    auto n = w.size();

    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto numEdges = vector<int>(n,0);
    for(int i=0; i<n; ++i)
        for(int j=i+1; j<n; ++j)
            if(distance(pos[i], pos[j]) < c*std::pow(w[i] * w[j] / W, 1.0/d)) {
                numEdges[i]++;
                numEdges[j]++;
            }

    auto res = std::accumulate(numEdges.begin(), numEdges.end(), 0.0);
    return res/n;
}


void test_estimation(int seed) {

    auto n = 1000;
    auto PLE = -2.5;
    auto desired_avg = 25;
    auto runs = 50;

    auto weights = generateWeights(n, PLE, seed);

    // do the tests for all dimensions < 5
    for(auto d = 1; d<5; ++d) {
        auto estimated_c = estimate_scaling(weights, desired_avg, d);

        // observed avg with estimated c (over multiple runs with different positions)
        auto observed_avg = 0.0;
        for(int i = 0; i<runs; ++i) {

            auto position_seed = i;

            // try GIRGS generator
            auto generator = Generator();
            if(estimated_c >= 1.0){
                auto scaled_weights = weights;
                auto scaling_factor = pow(estimated_c, d);
                for(auto& each : scaled_weights) // scale weights according to c
                    each = scaling_factor*each;
                generator.generateGIRG(d, scaled_weights, std::numeric_limits<double>::infinity(), 1.0, position_seed);
            } else {
                generator.generateGIRG(d, weights, std::numeric_limits<double>::infinity(), estimated_c, position_seed);
            }
            auto generated_avg_generator = generator.avg_degree();

            // try quadratic sampling with same positions as the generator
            auto pos = vector<vector<double>>(n);
            for(int i=0; i<n; ++i)
                pos[i] = generator.graph()[i].coord;
            auto generated_avg_trivial = quadraticSampling(weights,estimated_c, d, pos);

            // generator must yield same results as quadratic sampling
            test(generated_avg_generator == generated_avg_trivial);
            observed_avg += generated_avg_generator;
        }
        observed_avg /= runs;

        test(abs(desired_avg - observed_avg) < 0.1);
    }
}


int main(int argc, char* argv[]) {

    const auto seed = 1337;

    test_estimation(seed);


    cout << "all tests passed." << endl;
    return 0;
}
