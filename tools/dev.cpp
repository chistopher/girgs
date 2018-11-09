
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

#include <Generator.h>


using namespace std;


double trivialTryAvgWithPos(const std::vector<double>& w, double c, int d, const vector<vector<double>>& pos) {

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

vector<vector<double>> samplePositions(int n, int dimension, int positionSeed) {
    auto gen = std::mt19937(positionSeed >= 0 ?  positionSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)

    auto result = vector<vector<double>>(n);
    for(int i=0; i<n; ++i) {
        result[i].resize(dimension);
        for (int d=0; d<dimension; ++d)
            result[i][d] = dist(gen);
    }
    return result;
}


int main(int argc, char* argv[]) {

    auto d = 2;
    auto n = 1000;
    auto desired_avg = 9;
    auto alpha = numeric_limits<double>::infinity();

    auto weight_seed = 1337;
    auto PLE = -2.5;

    Generator generator;
    generator.setWeights(n, PLE, weight_seed);
    auto weights = generator.weights();

    auto start = std::chrono::high_resolution_clock::now();
    auto scaling = generator.scaleWeights(desired_avg, d, alpha);
    auto estimated_c = pow(scaling, 1.0/d);
    auto end = std::chrono::high_resolution_clock::now();

    int runs = 10;
    auto observed_avg = 0.0;
    auto time_sum = 0;
    for(int i = 0; i<runs; ++i) {

        auto positions = samplePositions(n, d, i);
        generator.setPositions(positions);

        auto start = std::chrono::high_resolution_clock::now();
        generator.generateTreshold();
        auto end = std::chrono::high_resolution_clock::now();
        time_sum += std::chrono::duration_cast<chrono::milliseconds>(end-start).count();

        auto avg_deg = generator.avg_degree();
        auto quadratic_avg = trivialTryAvgWithPos(weights, estimated_c, d, positions);
        if(avg_deg != quadratic_avg)
            cout << "oh no!" << endl;
        observed_avg += avg_deg;
    }
    observed_avg /= runs;

    cout << "find scaling (ms) " << std::chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
    cout << "generation (ms)   " << time_sum / runs << endl;
    cout << "desired avg deg   " << desired_avg << endl;
    cout << "observed avg deg  " << observed_avg << endl;
    
    return 0;
}
