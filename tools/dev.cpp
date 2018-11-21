
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

#include <Generator.h>


using namespace std;


double trivialTryAvgWithPos(const std::vector<double>& w, double c, int d, const vector<vector<double>>& pos, double a) {

    auto n = w.size();

    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto sum = 0.0;
    for(int i=0; i<n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (a == std::numeric_limits<double>::infinity()) {
                if (distance(pos[i], pos[j]) < c * std::pow(w[i] * w[j] / W, 1.0 / d)) {
                    sum += 2;
                }
            } else {
                sum += 2 * min(c / pow(distance(pos[i], pos[j]), a * d) * pow(w[i] * w[j] / W, a), 1.0);
            }
        }
    }

    return sum/n;
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
    auto desired_avg = 15;
    auto alpha = 4.0;

    auto weight_seed = 1337;
    auto position_base_seed = random_device()();
    auto PLE = -2.5;

    Generator generator;
    generator.setWeights(n, PLE, weight_seed);
    auto weights = generator.weights();

    auto scaling = generator.scaleWeights(desired_avg, d, alpha);
    auto estimated_c = alpha != std::numeric_limits<double>::infinity()
            ? pow(scaling, alpha)
            : pow(scaling, 1.0/d);


    int runs = 50;
    auto observed_avg = 0.0;
    auto expected_avg = 0.0;
    for(int i = 0; i<runs; ++i) {

        auto positions = samplePositions(n, d, position_base_seed+i);
        generator.setPositions(positions);

        generator.generate(alpha, i+n);

        auto s1 = chrono::high_resolution_clock::now();
        auto avg_deg = generator.avg_degree();
        auto s2 = chrono::high_resolution_clock::now();
        auto quadratic_avg = trivialTryAvgWithPos(weights, estimated_c, d, positions, alpha);
        auto e = chrono::high_resolution_clock::now();

        auto t1 = chrono::duration_cast<chrono::nanoseconds>(s2-s1).count();
        auto t2 = chrono::duration_cast<chrono::nanoseconds>(e-s2).count();
        cout  << "linear: " << t1 << "\tquadratic: " << t2 << endl;
        cout  << "observed: " << avg_deg << "\texprected: " << quadratic_avg << endl;
        observed_avg += avg_deg;
        expected_avg += quadratic_avg;
    }
    observed_avg /= runs;
    expected_avg /= runs;

    cout << endl << "SUMMARY: " << endl;
    cout << "desired average degree           " << desired_avg << endl;
    cout << "estimated constant               " << estimated_c << endl;
    cout << "equivalent weight scaling factor " << scaling << endl;

    cout << "expected average degree for estimated constant  " << expected_avg << endl;
    cout << "observed average degree using estimated constant " << observed_avg << endl;

    return 0;
}
