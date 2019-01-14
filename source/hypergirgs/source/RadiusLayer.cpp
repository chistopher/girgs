
#include <hypergirgs/RadiusLayer.h>

#include <cassert>

#include <hypergirgs/AngleHelper.h>


namespace hypergirgs {


RadiusLayer::RadiusLayer(double r_min, double r_max,
                         unsigned int targetLevel, unsigned int firstCell,
                         const Point* begin, const Point* end)
: m_r_min(r_min)
, m_r_max(r_max)
, m_target_level(targetLevel)
, m_base(begin)
{
    const auto cellsInLevel = AngleHelper::numCellsInLevel(targetLevel);

    m_prefix_sums.resize(cellsInLevel+1, 0);

    // count num of points in each cell
    int i = 0;
    for(auto it = begin; it != end; ++it, ++i) {
        // check radius lies within layer's radial bound
        assert(m_r_min <= it->radius);
        assert(it->radius <= m_r_max);

        // check cell id lies with layer's bound
        assert(firstCell <= it->cell_id);
        assert(it->cell_id < firstCell + cellsInLevel);

        const auto cell = it->cell_id - firstCell;
        ++m_prefix_sums[cell];
    }

    // compute exclusive prefix sums
    // prefix_sums[i] is the number of all points in cells j<i of the same level
    {
        assert(m_prefix_sums.back() == 0);
        unsigned sum = 0;
        for(auto& val : m_prefix_sums)  {
            const auto tmp = val;
            val = sum;
            sum += tmp;
        }
        assert(m_prefix_sums.back() == std::distance(begin, end));
    }
}

} // namespace hypergirgs
