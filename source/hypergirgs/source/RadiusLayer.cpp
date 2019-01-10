
#include <hypergirgs/RadiusLayer.h>

#include <cassert>

#include <hypergirgs/AngleHelper.h>


namespace hypergirgs {


RadiusLayer::RadiusLayer(double r_min, double r_max, unsigned int targetLevel, const std::vector<int> &nodes,
                         const std::vector<double> &angles, const std::vector<Point> &points)
: m_r_min(r_min)
, m_r_max(r_max)
, m_target_level(targetLevel)
{
    // allocate stuff
    const auto cellsInLevel = AngleHelper::numCellsInLevel(targetLevel);
    m_prefix_sums.resize(cellsInLevel+1, 0);

    // count num of points in each cell
    for(auto node : nodes) {
        const auto cell = AngleHelper::cellForPoint(angles[node], targetLevel);
        assert(cell < cellsInLevel);
        ++m_prefix_sums[cell];
    }

    // compute exclusive prefix sums
    // prefix_sums[i] is the number of all points in cells j<i of the same level
    {
        unsigned sum = 0;
        for(auto& val : m_prefix_sums)  {
            const auto tmp = val;
            val = sum;
            sum += tmp;
        }
    }

    m_points.resize(m_prefix_sums.back());

    // fill point lookup
    auto num_inserted = std::vector<int>(cellsInLevel, 0); // keeps track of bucket size for counting sort
    for(auto node : nodes){
        auto targetCell = AngleHelper::cellForPoint(angles[node], targetLevel);
        m_points[m_prefix_sums[targetCell] + num_inserted[targetCell]] = points[node];
        ++num_inserted[targetCell];
    }
}

} // namespace hypergirgs
