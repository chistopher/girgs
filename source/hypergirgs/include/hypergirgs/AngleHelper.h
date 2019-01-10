
#pragma once

#include <utility>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {


class HYPERGIRGS_API AngleHelper
{
public:

    // static helper functions
    static constexpr unsigned int numCellsInLevel(unsigned int level) noexcept { return 1u<<level; }
    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept { return (1u<<level)-1; }

    static constexpr unsigned int parent(unsigned int cell) noexcept { return (cell-1)/2; }
    static constexpr unsigned int firstChild(unsigned int cell) noexcept { return 2*cell+1; }
    static constexpr unsigned int secondChild(unsigned int cell) noexcept { return 2*(cell+1); }
    static constexpr unsigned int numChildren() noexcept { return 2; }

    static std::pair<double,double> bounds(unsigned int cell, unsigned int level);
    static unsigned int cellForPoint(double angle, unsigned int targetLevel);

    static bool touching(unsigned int cellA, unsigned int cellB, unsigned int level);

    // returns a lower bound for the angular difference of two points in these cells
    static double dist(unsigned int cellA, unsigned int cellB, unsigned int level);
};


} // namespace hypergirgs
