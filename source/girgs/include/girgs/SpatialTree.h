
#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <numeric>
#include <cassert>

#include <omp.h>

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
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);

    /**
     * @brief
     *  Same as visitCellPair(unsigned int, unsigned int, unsigned int) but stops recursion before first_parallel_level.
     *  Instead, the calls that would be made in this level are saved in parallel_calls.
     *  The saved calls are grouped by their (level local) cellA parameter.
     *
     * @param cellA
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param cellB
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param level
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param first_parallel_level
     *  The level before which we "saw off" the recursion.
     *  To get sufficient parallel cells (the outer size of parallel_calls) this should be computed as
     *  \f$ 2^{dl} \geq kt \f$ solved for l (d dimension, l first_parallel_level, t threads, k tuning parameter).
     *  We get \f$ l \geq \log_2(kt) / d \f$.
     * @param parallel_calls
     */
    void visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
            unsigned int first_parallel_level, std::vector<std::vector<unsigned int>>& parallel_calls);

    // sample edges of type 1 between V_i^A V_j^B
    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);
    // sample edges of type 2 between V_i^A V_j^B
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    // called v(i) in paper (=v(i,0)); not that our layers are offset by 1 (i in paper is i+1 for us)
    unsigned int weightLayerTargetLevel(int layer) const;
    // called v(i,j) in paper; not that our layers are offset by 1 (i in paper is i+1 for us)
    unsigned int partitioningBaseLevel(int layer1, int layer2) const;

    bool checkEdgeExplicit(double dist, double w1, double w2);

protected:

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
   
    std::vector<std::mt19937> m_gens; // for parallel access
    std::vector<std::uniform_real_distribution<>> m_dists; // for parallel access

#ifndef NDEBUG
    std::vector<long long> m_type1_checks; ///< number of node pairs per thread that are checked via a type 1 check
    std::vector<long long> m_type2_checks; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};


} // namespace girgs

#include <girgs/SpatialTree.inl>
