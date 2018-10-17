
#pragma once

#include <array>
#include <vector>
#include <utility>
#include <cassert>

// datastucture is full D-dimensional spatial tree with L levels
// TODO does each function need to be passed the level of the cell?
template<unsigned int D, unsigned int L>
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
    static const unsigned int numCells = firstCellOfLevel(L);
    static const unsigned int numChildren = 1u<<D;

    const static auto level = L;
    const static auto dimension = D;

    constexpr SpatialTree(); // constexpr c'tor

    // for testing purpose
    std::array<std::pair<double,double>, D> bounds(int cell, unsigned int level) const;
    // used during runtime for countingsort of points
    unsigned int cellForPoint(std::array<double, D>& point, unsigned int targetLevel) const;

    // graph generation stuff
    struct Node {
        std::array<double, D> coord;
        int weight;
        int index; // for debug // TODO maybe remove?
        std::vector<Node*> m_edges;
    };
    using Graph = std::vector<Node>;
    // entry point for the algorithm
    Graph generateGraph(std::vector<int>& weights);

protected:

    // recursive function that samples all edges between points in A and B
    void visitCellPair(unsigned int cellA, unsigned int cellB, Graph& graph) const;


    // helper functions from paper:

    // returns the number of points in a cell
    // PRECONDITION: cell must be of lower level than the target level
    //               i.e. cell index > first cell index in target level
    // begin = first descendent of cell in targetLevel
    // end = last child of cell in targetLevel
    // return how many points there are in cells {begin..end} using prefix sums
    int pointsInCell(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel) const; // maybe also give level of cell?

    // returns the k-th point in a cell
    // PRECONDITION: cell must be of lower level than the target level
    //               i.e. cell index > first cell index in target level (see n)
    // find first descendent of cell in target level and use m_A to find the point
    Node* kthPoint(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel, int k) const; // maybe also give level of cell?


private:
    // compile time data:

    std::array<std::array<int, D>, numCells> m_coords; // TODO flatten coords array like in inverse mapping
    std::array<int, numCells> m_coords2Index;

    // runtime data:

    // the number of points in each cell
    std::array<int, numCells> m_points_in_cell;
    // for each cell $i$ the number of points in all cells $j<i$ of the same level (prefix sums)
    std::array<int, numCells> m_prefix_sums;
    // m_A[l][m_prefix_sums[i]+k] contains the k-th point in the i-th cell of level l
    std::array<std::vector<Node*>, L>* m_A; // ptr to enable constexpr
    // breaks recursion after maxLevel so that a tree can be used for fewer layers than the tree has level
    unsigned int m_maxLevel;
    // minimum weigh in weight layer w_i
    std::vector<std::vector<Node*>>* m_weight_layers; // ptr to enable constexpr
};


#include <SpatialTree.inl>