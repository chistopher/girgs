
#pragma once

#include <array>
#include <vector>
#include <utility>
#include <cassert>


// datastucture is full D-dimensional spatial tree with L levels
template<unsigned int D, unsigned int L>
class SpatialTree
{
public:

    // the total number of cells on all levels {0..(L-1)}
    // $\sum_{i=0}^{L-1} 2^{DL} = \frac{2^{DL}-1}{2^D-1}
    static constexpr unsigned int numCellsInLevel(unsigned int level) { return 1<<(D*level); }
    static constexpr unsigned int firstCellOfLevel(unsigned int level) { return ((1<<(D*level))-1)/((1<<D)-1); }
    static constexpr unsigned int firstChild(unsigned int cell) { return (1<<D)*cell+1; }
    static constexpr unsigned int parent(unsigned int cell) { return (cell-1)/(1<<D); }
    static const unsigned int numCells = firstCellOfLevel(L);
    static const unsigned int numChildren = 1<<D;

    const static auto level = L;
    const static auto dimension = D;


    SpatialTree() {} // for now for testing
    // SpatialTree(std::array<std::vector<std::array<D>>, L>& weigh_levels); // TODO change int to pointer to actual points

    // returns the number of points in a cell
    // PRECONDITION: cell must be of lower level than the target level
    //               i.e. cell index > first cell index in target level (see n)
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
    // int level(int cell) const;
    std::array<int, D> coords(int cell) const;
    std::array<std::pair<double,double>, D> bounds(int cell) const;

    int parent(int cell) const;
    int firstChild(int cell) const;

    int cellByIndex(std::array<int, D>, int level) const;
    int cellByCoord(std::array<double, D>, int level) const;

private:
    // the number of points in each cell
    std::array<int, numCells> m_points_in_cell;
    // for each cell $i$ the number of points in all cells $j<i$ of the same level (prefix sums)
    std::array<int, numCells> m_prefix_sums;

    // m_A[l][m_prefix_sums[i]+k] contains the k-th point in the i-th cell of level l
    std::array<std::vector<int>, L> m_A;
};


template<unsigned int D, unsigned int L>
int SpatialTree<D,L>::pointsInCell(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel) const {
    assert(fromLevel <= targetLevel);
    assert(firstCellOfLevel(fromLevel) <= cell && cell < firstCellOfLevel(fromLevel+1)); // cell is from fromLevel

    // we want the begin-th and end-th cell in level targetLevel to be the first and last descendent of cell in this level
    auto begin = 0;
    auto end = 0;

    // first descendent
    /* OP1: O(L)
    begin = cell;
    for(int i=fromLevel; fromLevel<targetLevel; fromLevel++)
        begin = firstChild(begin);
    */
    // OPT2: O(1) (use local index in layer instead of global index)
    auto localIndexCell = cell - firstCellOfLevel(fromLevel);
    auto localIndexDescendent = localIndexCell * (1<<D)<<(targetLevel - fromLevel); // each cell before the parent splits in 2^D cells in the next layer that are all before our descendent
    begin = localIndexDescendent + firstCellOfLevel(targetLevel);

    // last descendent
    // do OPT2 to determine first descendent of next cell and subtract one to get our last descendent
    auto localIndexCellNext = localIndexCell + 1;
    auto localIndexDescendentOfNext = (localIndexCellNext) * (1<<D)<<(targetLevel - fromLevel);
    end = localIndexDescendentOfNext - 1 + firstCellOfLevel(targetLevel);

    assert(begin < firstCellOfLevel(targetLevel+1));
    assert(end < firstCellOfLevel(targetLevel+1));

    std::cout << cell << '\t' << fromLevel << std::endl;
    std::cout << begin << '\t' << end << '\t' << targetLevel << std::endl;
    return m_prefix_sums[end] - m_prefix_sums[begin] + m_points_in_cell[end];
}
