
#pragma once

#include <vector>


namespace hypergirgs {


class AngleHelper
{
public:

    // static helper functions
    static constexpr unsigned int numCellsInLevel(unsigned int level) noexcept { return 1u<<level; }
    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept { return (1u<<level)-1; }

    static constexpr unsigned int parent(unsigned int cell) noexcept { return (cell-1)/2; }
    static constexpr unsigned int firstChild(unsigned int cell) noexcept { return 2*cell+1; }
    static constexpr unsigned int secondChild(unsigned int cell) noexcept { return 2*(cell+1); }
    static constexpr unsigned int numChildren() noexcept { return 2; }

    AngleHelper() = default;
    explicit AngleHelper(unsigned int levels);

    std::pair<double,double> bounds(unsigned int cell, unsigned int level) const;
    unsigned int cellForPoint(double angle, unsigned int targetLevel) const;

    bool touching(unsigned int cellA, unsigned int cellB, unsigned int level) const;

    // returns a lower bound for the angular difference of two points in these cells
    double dist(unsigned int cellA, unsigned int cellB, unsigned int level) const;

    unsigned int levels() const { return m_levels; }


protected:
    unsigned int m_levels = 0;

    std::vector<int> m_coords;       // cell index to coord
    std::vector<unsigned int>       m_coords2Index; // packed coord to cell index
};


} // namespace hypergirgs
