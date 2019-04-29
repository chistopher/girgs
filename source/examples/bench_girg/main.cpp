
#include <iostream>
#include <cassert>
#include <limits>

#include <omp.h>

#include <girgs/Generator.h>
#include <girgs/SpatialTree.h>
#include <girgs/ScopedTimer.h>

#include "CounterPerThread.h"


using namespace std;
using namespace girgs;

// returns number of edges
template<unsigned int D>
void measure_edge_generation(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, double alpha, int samplingSeed,
        CounterPerThread<uint64_t>& counter, double& time_pre, double& time_edges) {
    // Preprocess
    auto addEdge = [&counter] (int, int, int tid) {counter.add(tid);};
    auto generator = [&] {
        ScopedTimer timer("", time_pre);
        return makeSpatialTree<D>(weights, positions, alpha, addEdge);
    }();
    // generate
    {
        ScopedTimer timer("", time_edges);
        generator.generateEdges(samplingSeed);
    }
}


void measure(int dimension, int n, int avgDeg, double alpha, double ple, int threads, int seed, int plot) {

    omp_set_num_threads(threads);
    assert(threads == omp_get_max_threads());

    auto weightSeed = seed + 1000;
    auto positionSeed = seed+ 10000;
    auto samplingSeed = seed+ 100000;

    double time_weights, time_positions, time_binary, time_pre, time_edges, time_total;

    CounterPerThread<uint64_t> counter_num_edges;
    {
        ScopedTimer total_timer("", time_total);

        auto weights = [&]{
            ScopedTimer timer("", time_weights);
            return generateWeights(n, ple, weightSeed);
        }();

        auto positions = [&] {
            ScopedTimer timer("", time_positions);
            return generatePositions(n, dimension, positionSeed);
        }();

        {
            ScopedTimer timer("", time_binary);
            scaleWeights(weights, avgDeg, dimension, alpha);
        }

        {
            switch(dimension) {
                case 1: measure_edge_generation<1>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
                case 2: measure_edge_generation<2>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
                case 3: measure_edge_generation<3>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
                case 4: measure_edge_generation<4>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
                case 5: measure_edge_generation<5>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
                default: measure_edge_generation<1>(weights, positions, alpha, samplingSeed, counter_num_edges, time_pre, time_edges); break;
            }
        }
    }

    auto edges = counter_num_edges.total();
    auto degree = 2.0 * edges / n;

    cout << dimension << ','
         << n << ','
         << avgDeg << ','
         << alpha << ','
         << ple << ','
         << threads << ','
         << seed << ','
         << plot << ','
         << time_weights << ','
         << time_positions << ','
         << time_binary << ','
         << time_pre << ','
         << time_edges << ','
         << time_total << ','
         << edges << ','
         << degree << '\n';
}


int main(int argc, char* argv[]) {

    cout << "dimension,n,avgDeg,alpha,ple,threads,seed,plot,TimeWeights,TimePositions,TimeBinary,TimePre,TimeEdges,TimeTotal,GenNumEdge,GenAvgDeg\n";

    int seed = 0;

    auto d = 1;
    auto n = 1<<15;
    auto alpha = std::numeric_limits<double>::infinity();
    auto alpha_binomial = 2.0;
    auto ple = 2.5;
    auto deg=10;
    auto threads = 1;
    auto reps = 10;

    for(int rep=0; rep<reps; ++rep) {
        clog << "rep " << rep << endl;

        clog << "growing n" << endl;
        for(int i = 1<<10; i<= (1<<20); i <<= 1) {
            clog << i << endl;
            for(auto deg : {10,100,1000})
		if(deg*4 < i) 
		    measure(d, i, deg, alpha, ple, threads, ++seed, 0);
        }


        clog << "growing deg" << endl;
        for(int i = 2; i<= 1<<12; i <<= 1) {
            clog << i << endl;
            measure(d, n, i, alpha, ple, threads, ++seed, 1);
        }

        clog << "shrinking alpha" << endl;
        for(auto i : {alpha, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0}) {
            clog << i << endl;
            measure(d, n, deg, i, ple, threads, ++seed, 2);
        }

        clog << "growing dimension" << endl;
        for(auto i : {1,2,3,4,5}) {
            clog << i << endl;
            measure(i, n, deg, alpha, ple, threads, ++seed, 3);
        }
        /*
        clog << "strong scaling threshold" << endl;
        for(auto i : {1,2,3,4,5,6}) {
            clog << i << endl;
            for(auto d : {1,2,3,4,5}) {
                measure(d, n, deg, alpha, ple, i, ++seed, 4);
            }
        }

        clog << "strong scaling binomial" << endl;
        for(auto i : {1,2,3,4,5,6}) {
            clog << i << endl;
            for(auto d : {1,2,3,4,5}) {
                measure(d, n, deg, alpha_binomial, ple, i, ++seed, 5);
            }
        }
        */
    }

    return 0;
}
