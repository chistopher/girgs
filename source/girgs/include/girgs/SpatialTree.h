
#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <numeric>
#include <cassert>

#include <girgs/SpatialTreeCoordinateHelper.h>
#include <girgs/WeightLayer.h>
#include <girgs/Node.h>


namespace girgs {


/**
 * @brief
 *  the data stucture is a full D-dimensional space partitioning tree
 */
template<unsigned int D>
class SpatialTree
{
public:
    static const auto dimension = D;

    SpatialTree() = default;

    // entry point for the algorithm
    void generateEdges(std::vector<Node>& graph, double alpha, int seed);

protected:

    // recursive function that samples all edges between points in A and B
    void visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
            unsigned int first_parallel_level, std::vector<std::vector<unsigned int>>& parallel_calls);
    void visitRoot_parallel();
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);

    // sample edges of type 1 between V_i^A V_j^B
    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);
    // sample edges of type 2 between V_i^A V_j^B
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    // called v(i) in paper
    unsigned int weightLayerTargetLevel(int layer) const;
    // called v(i,j) in paper
    unsigned int partitioningBaseLevel(int layer1, int layer2) const;

    bool checkEdgeExplicit(double dist, double w1, double w2);

private:

    unsigned int m_layers; // number of layers
    unsigned int m_levels; // number of levels

    SpatialTreeCoordinateHelper<D> m_helper;        // computes index to coordinate mappings, also offers static helper for cell indices
    std::vector<WeightLayer<D>> m_weight_layers;    // stores all nodes of one weight layer and provides the data structure described in paper
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs; // which pairs of weight layers to check in each level

    double m_w0;                // minimum weight
    double m_wn;                // maximum weight
    double m_W;                 // sum of weights
    int    m_baseLevelConstant; // log(W/w0^2)

    double m_alpha;             // edge prob stuff
    std::mt19937 m_gen;
    std::uniform_real_distribution<> m_dist;

    long int m_1, m_2;          // number of checked type 1 and type 2 connections
};


} // namespace girgs

#include <girgs/SpatialTree.inl>
