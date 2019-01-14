
#include <chrono>
#include <iostream>

#include <girgs/Generator.h>


using namespace std;
using namespace girgs;

// GIRG parameter
const auto ple = -2.5;
const auto alpha = 1.1;
//const auto alpha = std::numeric_limits<double>::infinity();
const auto avgDeg = 10;

// measuring parameter
const auto SEED = 13;
const auto RUNS = 10;
const auto N_1 = 12500  * (alpha == std::numeric_limits<double>::infinity() ? 10 : 1);
const auto N_2 = 100100 * (alpha == std::numeric_limits<double>::infinity() ? 10 : 1);
const auto D_MAX = 4;
using resolution = chrono::milliseconds;


int measure(int dimension, int n, int run, Generator& g) {

    auto positionSeed = SEED+run;
    auto samplingSeed = SEED+run+n;

    g.setPositions(n, dimension, positionSeed);
    
    // measure
    auto start = chrono::high_resolution_clock::now();
    g.generate(alpha, samplingSeed);
    auto end = chrono::high_resolution_clock::now();

    return chrono::duration_cast<resolution>(end-start).count();
}

int avg(int dimension, int n){

	// do scaling outside of measurement loop
	Generator g;
	g.setWeights(n, ple, SEED);
	g.scaleWeights(avgDeg, dimension, alpha);

    auto sum = 0;
    for (int i = 0; i < RUNS; ++i)
        sum = sum + measure(dimension, n, i, g);
    return sum / RUNS;
}


int main(int argc, char* argv[]) {

    for(auto d = 1; d<D_MAX; ++d){
        cout << "d=" << d << '\n';
        for(auto n = N_1; n<N_2; n*=2){
            auto ggg =  avg(d,n);
            cout << n << '\t' << ggg << endl;
        }
        cout << endl;
    }
    return 0;
}
