
#pragma once

#include <vector>
#include <cmath>

#include <hypergirgs/AngleHelper.h>


namespace hypergirgs {


class RadiusLayer {
public:

	RadiusLayer() = delete;

	RadiusLayer(unsigned int layer, unsigned int targetLevel, const AngleHelper &helper,
				std::vector<int> &&nodes, const std::vector<double> &radii);


    int pointsInCell(unsigned int cell, unsigned int level) const;

    int kthPoint(unsigned int cell, unsigned int level, int k) const;

    const int * firstPointPointer(unsigned int cell, unsigned int level) const;

protected:

    const unsigned int m_layer;             ///< the index of the layer
    const unsigned int m_target_level;      ///< the insertion level for the current weight layer (v(i) = wiw0/W)

    std::vector<int>   m_points_in_cell;    ///< the number of points in each cell of target_level
    std::vector<int>   m_prefix_sums;       ///< for each cell c in target level: the sum of points of this layer in all cells <c
    std::vector<int>   m_A;                 ///< m_A[m_prefix_sums[i]+k] contains the k-th point in the i-th cell of target level
};



} // namespace hypergirgs
