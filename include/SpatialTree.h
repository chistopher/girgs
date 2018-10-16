
#pragma once

#include <array>
#include <vector>
#include <utility>
#include <optional>
#include <cassert>


// datastucture is full D-dimensional spatial tree with L levels
// TODO does each function need to be passed the level of the cell?
template<unsigned int D, unsigned int L>
class SpatialTree
{
public:

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


    constexpr SpatialTree(); // for now for testing
    // SpatialTree(std::array<std::vector<std::array<D>>, L>& weigh_levels); // TODO change int to pointer to actual points

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
    int kthPoint(unsigned int cell, int k, unsigned int targetLevel) const; // maybe also give level of cell?


    // for testing purpose
    // maybe need them later
    std::array<std::pair<double,double>, D> bounds(int cell, unsigned int level) const;
    unsigned int cellForPoint(std::array<double, D>& point, unsigned int targetLevel) const;

private:
    std::array<std::array<int, D>, numCells> m_coords; // TODO flatten coords array like in inverse mapping
    std::array<int, numCells> m_coords2Index;

    // the number of points in each cell
    std::array<int, numCells> m_points_in_cell;
    // for each cell $i$ the number of points in all cells $j<i$ of the same level (prefix sums)
    std::array<int, numCells> m_prefix_sums;
    // m_A[l][m_prefix_sums[i]+k] contains the k-th point in the i-th cell of level l
    std::array<std::vector<int>, L>* m_A;
};


#include <SpatialTree.inl>