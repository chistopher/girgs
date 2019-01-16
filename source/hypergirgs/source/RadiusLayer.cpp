
#include <hypergirgs/RadiusLayer.h>

#include <cassert>

#include <hypergirgs/AngleHelper.h>


namespace hypergirgs {

RadiusLayer::RadiusLayer(double r_min, double r_max, unsigned int targetLevel,
            const Point* base, const unsigned int* prefix_sum)
    : m_r_min{r_min}, m_r_max{r_max}, m_target_level{targetLevel},
      m_base{base}, m_prefix_sums{prefix_sum}
{
#ifndef NDEBUG
    const auto cellsInLevel = AngleHelper::numCellsInLevel(targetLevel);

    const auto begin = m_base + m_prefix_sums[0];
    const auto end   = m_base + m_prefix_sums[cellsInLevel];

    for(auto it = begin; it != end; ++it) {
        // check radius lies within layer's radial bound
        assert(m_r_min <= it->radius);
        assert(it->radius <= m_r_max);
    }
#endif
}
} // namespace hypergirgs
