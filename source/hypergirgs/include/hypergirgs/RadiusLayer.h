
#pragma once

#include <cassert>
#include <vector>
#include <utility>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/Point.h>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {


class HYPERGIRGS_API RadiusLayer {
public:

	RadiusLayer() = delete;

	RadiusLayer(double r_min, double r_max, unsigned int targetLevel,
                const Point* base,
                const unsigned int* prefix_sum);

    int pointsInCell(unsigned int cell, unsigned int level) const {
        auto cellBoundaries = levelledCell(cell, level);
        assert(cellBoundaries.first  + AngleHelper::firstCellOfLevel(level) < AngleHelper::firstCellOfLevel(m_target_level+1));
        assert(cellBoundaries.second + AngleHelper::firstCellOfLevel(level) < AngleHelper::firstCellOfLevel(m_target_level+1));

        return m_prefix_sums[cellBoundaries.second+1] - m_prefix_sums[cellBoundaries.first];
    }

    const Point& kthPoint(unsigned int cell, unsigned int level, int k) const {
        auto cellBoundaries = levelledCell(cell, level);
        return m_base[m_prefix_sums[cellBoundaries.first] + k];
    }

    std::pair<const Point*, const Point*> cellIterators(unsigned int cell, unsigned int level) const {
        auto cellBoundaries = levelledCell(cell, level);
        const auto begin_end = std::make_pair(
                m_base + m_prefix_sums[cellBoundaries.first],
                m_base + m_prefix_sums[cellBoundaries.second+1]);
        assert(begin_end.first <= begin_end.second);
        return begin_end;
    }

    // static generation and helper
    static std::vector<RadiusLayer>
    buildPartition(const std::vector<double>& radii, const std::vector<double>& angles,
                   const double R, const double layer_height,
                   std::vector<Point>& points, std::vector<unsigned int>& first_in_cell, // output parameter
                   bool enable_profiling);


    // takes lower bound on radius for two layers
    static unsigned int partitioningBaseLevel(double r1, double r2, double R) noexcept {
        assert(r1 < R && r2 < R);

        auto level = 0u;
        auto cellDiameter = 2.0*PI;
        // find deepest level in which points in all non-touching cells are not connected
        while(hypergirgs::hyperbolicDistance(r1, 0, r2, (cellDiameter/2)) > R){
            level++;
            cellDiameter /= 2;
        }
        return level;
    }

protected:

    std::pair<unsigned int, unsigned int> levelledCell(unsigned int cell, unsigned int level) const {
        assert(level <= m_target_level);
        assert(AngleHelper::firstCellOfLevel(level) <= cell && cell < AngleHelper::firstCellOfLevel(level + 1)); // cell is from fromLevel

        // we want the begin-th and end-th cell in level targetLevel to be the first and last descendant of cell in this level
        // we could apply the firstChild function to find the first descendant but this is in O(1)
        auto descendants = AngleHelper::numCellsInLevel(m_target_level - level);
        auto localIndexCell = cell - AngleHelper::firstCellOfLevel(level);
        auto localIndexDescendant = localIndexCell * descendants; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendant
        auto begin = localIndexDescendant;
        auto end = begin + descendants - 1;

        assert(begin <= end);

        return {begin, end};
    }

public:
    const double m_r_min;               ///< lower bound on radius for points in this layer
    const double m_r_max;               ///< upper bound on radius for points in this layer
    const unsigned int m_target_level;  ///< insertion level for this radius layer

protected:
    const Point* m_base;                ///< sorted array of all points
    const unsigned int* m_prefix_sums;  ///< for each cell c in target level: sum of points in m_base before first node in c

};

} // namespace hypergirgs
