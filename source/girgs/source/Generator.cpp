#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <mutex>

#include <omp.h>

#include <girgs/Generator.h>
#include <girgs/SpatialTree.h>
#include <girgs/WeightScaling.h>


namespace girgs {

std::vector<double> generateWeights(int n, double ple, int weightSeed) {
    auto result = std::vector<double>(n);
    auto gen = std::mt19937(weightSeed >= 0 ?  weightSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        result[i] = std::pow((std::pow(n,-ple+1)-1)*dist(gen) + 1, 1/(-ple+1));
    return result;
}

std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed) {
    auto result = std::vector<std::vector<double>>(n, std::vector<double>(dimension));
    auto gen = std::mt19937(positionSeed >= 0 ?  positionSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        for (int d=0; d<dimension; ++d)
            result[i][d] = dist(gen);
    return result;
}

double scaleWeights(std::vector<double>& weights, double desiredAvgDegree, int dimension, double alpha) {
    // estimate scaling with binary search
    double scaling;
    if(alpha > 10.0)
        scaling = estimateWeightScalingThreshold(weights, desiredAvgDegree, dimension);
    else if(alpha > 0.0 && alpha != 1.0)
        scaling = estimateWeightScaling(weights, desiredAvgDegree, dimension, alpha);
    else
        throw("I do not know how to scale weights for desired alpha :(");

    // scale weights
    for(auto& each : weights) each *= scaling;
    return scaling;
}

std::vector<std::pair<int, int>> generateEdges(const std::vector<double> &weights, const std::vector<std::vector<double>> &positions,
        double alpha, int samplingSeed) {

    using edge_vector = std::vector<std::pair<int, int>>;
    edge_vector result;

    std::vector<std::pair<
            edge_vector,
            uint64_t[31] /* avoid false sharing */
    > > local_edges(omp_get_max_threads());

    constexpr auto block_size = size_t{1} << 20;

    std::mutex m;
    auto flush = [&] (const edge_vector& local) {
        std::lock_guard<std::mutex> lock(m);
        result.insert(result.end(), local.cbegin(), local.cend());
    };

    auto addEdge = [&](int u, int v, int tid) {
        auto& local = local_edges[tid].first;
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
        }
    };

    auto dimension = positions.front().size();

    switch(dimension) {
        case 1: makeSpatialTree<1>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 2: makeSpatialTree<2>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 3: makeSpatialTree<3>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 4: makeSpatialTree<4>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 5: makeSpatialTree<5>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        default:
            std::cout << "Dimension " << dimension << " not supported." << std::endl;
            std::cout << "No edges generated." << std::endl;
            break;
    }

    for(const auto& v : local_edges)
        flush(v.first);

    return result;
}


void saveDot(const std::vector<double> &weights, const std::vector<std::vector<double>> &positions,
             std::vector<std::pair<int, int>> graph, std::string file) {

    auto f = std::ofstream(file);
    f << "graph girg {\n\toverlap=scale;\n\n";
    for (int i = 0; i < weights.size(); ++i) {
        f << '\t' << i << " [label=\""
          << std::setprecision(2) << std::fixed << weights[i] << std::defaultfloat << std::setprecision(6)
          << "\", pos=\"";
        for (auto d = 0u; d < positions[i].size(); ++d)
            f << (d == 0 ? "" : ",") << positions[i][d];
        f << "\"];\n";
    }
    f << '\n';
    for (auto &edge : graph)
        f << '\t' << edge.first << "\t-- " << edge.second << ";\n";
    f << "}\n";
}

} // namespace girgs
