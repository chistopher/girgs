
#pragma once

#include <vector>
#include <array>
#include <cassert>


template<unsigned int D>
class SpatialTreeCoordinateHelper
{
public:

    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept { return ((1u<<(D*level))-1)/((1<<D)-1); }
    static constexpr unsigned int parent(unsigned int cell) noexcept { return (cell-1)/(1<<D); }
    static const auto numChildren = 1u<<D;

    // TODO move constructor
    SpatialTreeCoordinateHelper() = default;
    explicit SpatialTreeCoordinateHelper(unsigned int levels);

    std::array<std::pair<double,double>, D> bounds(unsigned int cell, unsigned int level) const;
    unsigned int cellForPoint(std::array<double, D>& point, unsigned int targetLevel) const;

    bool touching(unsigned int cellA, unsigned int cellB, unsigned int level) const;


    unsigned int levels() const { return m_levels; }

private:

    unsigned int m_levels = 0;

    std::vector<std::array<int, D>> m_coords; // cell index to coord TODO flatten coords array like in inverse mapping
    std::vector<unsigned int>       m_coords2Index; // packed coord to cell index
};


#include <SpatialTreeCoordinateHelper.inl>