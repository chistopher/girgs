
#pragma once

#include <utility>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {

constexpr double PI = 3.14159265358979323846;

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

    //! returns level local index of cell at angle
    static unsigned int cellForPoint(double angle, unsigned int targetLevel) {
        return static_cast<unsigned int>(angle / 2 / PI * numCellsInLevel(targetLevel));
    }

    static bool touching(unsigned int cellA, unsigned int cellB, unsigned int level);

    // returns a lower bound for the angular difference of two points in these cells
    static double dist(unsigned int cellA, unsigned int cellB, unsigned int level);

    static int cellsBetween(unsigned int cellA, unsigned int cellB, unsigned int level);
};


} // namespace hypergirgs
