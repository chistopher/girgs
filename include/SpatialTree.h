
#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <iostream>
#include <numeric>
#include <cassert>

#include <SpatialTreeCoordinateHelper.h>
#include <WeightLayer.h>


// datastucture is full D-dimensional space partitioning tree
template<unsigned int D>
class SpatialTree
{
public:

    // static helper functions
    // the total number of cells on all levels {0..(L-1)}
    // $\sum_{i=0}^{L-1} 2^{DL} = \frac{2^{DL}-1}{2^D-1}
    static constexpr unsigned int numCellsInLevel(unsigned int level) noexcept { return 1u<<(D*level); }
    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept { return ((1u<<(D*level))-1)/((1<<D)-1); }
    static constexpr unsigned int firstChild(unsigned int cell) noexcept { return (1<<D)*cell+1; }
    static constexpr unsigned int parent(unsigned int cell) noexcept { return (cell-1)/(1<<D); }

    static const auto numChildren = 1u<<D;
    static const auto dimension = D;

    SpatialTree() = default;

    using Graph = std::vector<Node>;
    // entry point for the algorithm
    void generateEdges(std::vector<Node>& graph, double alpha, int seed);

protected:

    // recursive function that samples all edges between points in A and B
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, Graph& graph);

    // sample edges of type 1 between V_i^A V_j^B
    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, Graph& graph);
    // sample edges of type 2 between V_i^A V_j^B
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, Graph& graph);

    // called v(i) in paper
    unsigned int weightLayerTargetLevel(int layer) const;
    // called v(i,j) in paper
    unsigned int partitioningBaseLevel(int layer1, int layer2) const;

    bool checkEdgeExplicit(double dist, double w1, double w2);

private:

    unsigned int m_layers = -1; // number of layers
    unsigned int m_levels = -1; // number of levels

    // can be precomputed at compile time
    SpatialTreeCoordinateHelper<D> m_helper;

    std::vector<WeightLayer<D>> m_weight_layers;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs; // which pairs of weight layers to check in each level

    double m_w0; // minimum weight
    double m_wn; // maximum weight
    double m_W; // sum of weights
    int    m_baseLevelConstant; // log(W/w0^2)

    // edge prob stuff
    double m_alpha;

    std::mt19937 m_gen;
    std::uniform_real_distribution<> m_dist;

    long int m_1, m_2; // test
};


#include <SpatialTree.inl>