
#pragma once

#include <vector>
#include <array>
#include <cassert>


namespace girgs {


template<unsigned int D>
class SpatialTreeCoordinateHelper
{
public:

    // static helper functions
    // the total number of cells on all levels {0..(L-1)}
    // $\sum_{i=0}^{L-1} 2^{DL} = \frac{2^{DL}-1}{2^D-1}
    static constexpr unsigned int numCellsInLevel(unsigned int level) noexcept { return 1u<<(D*level); }
    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept { return ((1u<<(D*level))-1)/((1<<D)-1); }

    static constexpr unsigned int parent(unsigned int cell) noexcept { return (cell-1)/(1<<D); }
    static constexpr unsigned int firstChild(unsigned int cell) noexcept { return (1<<D)*cell+1; }
    static constexpr unsigned int lastChild(unsigned int cell) noexcept { return firstChild(cell) + numChildren() - 1; }
    static constexpr unsigned int numChildren() noexcept { return 1u<<D; }

    SpatialTreeCoordinateHelper() = default;
    explicit SpatialTreeCoordinateHelper(unsigned int levels);

    std::array<std::pair<double,double>, D> bounds(unsigned int cell, unsigned int level) const;
    unsigned int cellForPoint(const std::vector<double>& point, unsigned int targetLevel) const;

    bool touching(unsigned int cellA, unsigned int cellB, unsigned int level) const;

    // implements the chebyshev distance metric (L_\infty)
    static double dist(std::vector<double>& a, std::vector<double>& b);

    // returns a lower bound for the distance of two points in these cells
    double dist(unsigned int cellA, unsigned int cellB, unsigned int level) const;


    unsigned int levels() const { return m_levels; }


protected:
    unsigned int m_levels = 0;

    std::vector<std::array<int, D>> m_coords;       // cell index to coord
};


} // namespace girgs

#include <girgs/SpatialTreeCoordinateHelper.inl>
