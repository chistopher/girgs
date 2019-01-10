
#pragma once

#include <vector>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {


class HYPERGIRGS_API RadiusLayer {
public:

	RadiusLayer() = delete;

	RadiusLayer(double r_min, double r_max, unsigned int targetLevel,
				const std::vector<int> &nodes, const std::vector<double> &angles);


    int pointsInCell(unsigned int cell, unsigned int level) const;

    int kthPoint(unsigned int cell, unsigned int level, int k) const;

    const int * firstPointPointer(unsigned int cell, unsigned int level) const;

public:
    const double m_r_min;
    const double m_r_max;
    const unsigned int m_target_level;

protected:
    std::vector<int> m_points_in_cell;    ///< the number of points in each cell of target_level
    std::vector<int> m_prefix_sums;       ///< for each cell c in target level: the sum of points of this layer in all cells <c
    std::vector<int> m_A;                 ///< m_A[m_prefix_sums[i]+k] contains the k-th point in the i-th cell of target level
};



} // namespace hypergirgs
