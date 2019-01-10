
#pragma once

#include <vector>

#include <hypergirgs/Point.h>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {


class HYPERGIRGS_API RadiusLayer {
public:

	RadiusLayer() = delete;

	RadiusLayer(double r_min, double r_max, unsigned int targetLevel,
				const std::vector<int> &nodes, const std::vector<double> &angles,
				const std::vector<Point> &points);


    int pointsInCell(unsigned int cell, unsigned int level) const;

    const Point& kthPoint(unsigned int cell, unsigned int level, int k) const;

    const Point* firstPointPointer(unsigned int cell, unsigned int level) const;

public:
    const double m_r_min;
    const double m_r_max;
    const unsigned int m_target_level;

protected:
    std::vector<int>   m_points_in_cell;    ///< the number of points in each cell of target_level
    std::vector<int>   m_prefix_sums;       ///< for each cell c in target level: the sum of points of this layer in all cells <c
    std::vector<Point> m_points; 			///< vector of points in this layer
};

} // namespace hypergirgs
