
#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>

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


double avg_degree(int n, double c, int dimension, double W, double sq_W, vector<double>& sorted_rich_club) {
    // compute overestimation
    auto overestimation = pow(2, dimension) * pow(c, dimension) / n * (W - sq_W/W);

    // subtract error
    auto error = 0.0;
    for(int i = 0; i<sorted_rich_club.size(); ++i)
        for(int j = 0; j<sorted_rich_club.size(); ++j) {
            if(i==j) continue;
            auto w1 = sorted_rich_club[i];
            auto w2 = sorted_rich_club[j];
            auto e = max( std::pow(2*c,dimension)*(w1*w2/W)-1.0, 0.0);
            error += e;
            if(e <= 0) break;
        }

    return overestimation - error/n;
}

double c_for_avg_deg(const std::vector<double>& weights, int desired_avg_degree, int dimension) {

    // compute some constant stuff
    auto max_weight = *max_element(weights.begin(), weights.end());
    auto W = 0.0, sq_W = 0.0;
    for(auto each : weights){
        W += each;
        sq_W += each*each;
    }

    // my function to do the exponential search on
    auto f = [W, sq_W, &weights, dimension, n(weights.size()), max_weight](double c) {
        vector<double> rich_club;
        for(auto weight : weights)
            if(std::pow(2*c,dimension) * (weight*max_weight/W) > 1.0)
                rich_club.push_back(weight);
        sort(rich_club.begin(), rich_club.end(), greater<>());
        return avg_degree(n, c, dimension, W, sq_W, rich_club);
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



int main(int argc, char* argv[]) {

    auto d = 1;
    auto n = 1000;
    auto desired_avg = 12;
    auto alpha = std::numeric_limits<double>::infinity();

    auto weight_seed = 1337;
    auto PLE = -2.5;
    auto weights = generateWeights(n, PLE, weight_seed);
    auto estimated_c = c_for_avg_deg(weights, desired_avg, d);

    int runs = 50;

    auto observed_avg = 0.0;

    for(int i = 0; i<runs; ++i) {

        auto position_seed = i;

        auto scaled_weights = weights;
        for(auto& each : scaled_weights) each *= estimated_c;

        auto generator = Generator();
        generator.generateGIRG(d, scaled_weights, alpha, 1.0, position_seed);

        auto pos = vector<vector<double>>(n);
        for(int i=0; i<n; ++i)
            pos[i] = generator.graph()[i].coord;

        auto generated_avg1 = generator.avg_degree();
        auto generated_avg2 =trivialTryAvgWithPos(weights,estimated_c, d, pos);
        if(generated_avg1 != generated_avg2)
            cout << "mismatch! " << generated_avg1 << " " << generated_avg2 << endl;
        observed_avg += generated_avg1;
    }
    observed_avg /= runs;


    cout << "desired avg deg  " << desired_avg << endl;
    cout << "estimated c      " << estimated_c << endl;
    cout << "observed avg deg " << observed_avg << endl;
    
    return 0;
}
