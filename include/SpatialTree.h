
#pragma once

#include <vector>
#include <algorithm>
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

    using Graph = std::vector<Node<D>>;
    // entry point for the algorithm
    Graph generateGraph(std::vector<double>& weights);

protected:

    // recursive function that samples all edges between points in A and B
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, Graph& graph) const;

    void sampleBetweenViAVjB(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, Graph& graph) const;


private:

    unsigned int m_layers; // number of layers
    unsigned int m_levels; // number of levels

    // can be precomputed at compile time
    SpatialTreeCoordinateHelper<D> m_helper;

    std::vector<WeightLayer<D>> m_weight_layers;

    double m_W; // sum of weights
};


#include <SpatialTree.inl>