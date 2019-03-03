#include <girgs/BitInterleaving.h>

namespace girgs {


template<unsigned int D>
SpatialTreeCoordinateHelper<D>::SpatialTreeCoordinateHelper(unsigned int levels)
    : m_levels(levels)
    , m_coords(firstCellOfLevel(levels))
{

    // calculate coords // TODO find implicit way to get from index to coords and back

    m_coords[0].fill(0); // first cell has [0]^D
    for(auto l=1u; l<levels; ++l) {
        for(auto cell = firstCellOfLevel(l); cell < firstCellOfLevel(l+1); ++cell) {
            m_coords[cell] = m_coords[parent(cell)];
            auto bitmask = ((cell-1)&(numChildren()-1)); // & is unnecessary but removes all but the last D bits
            for(auto d=0u; d<D; ++d) { // check orientation to parent in all dimensions with d-th bit in bitmask
                m_coords[cell][d] *= 2; // every level doubles each dimension
                m_coords[cell][d] += static_cast<bool>(bitmask&(1<<d)); // add 1 if we are "right" of parent in dimension d
            }
        }
    }
}


template<unsigned int D>
std::array<std::pair<double, double>, D> SpatialTreeCoordinateHelper<D>::bounds(unsigned int cell, unsigned int level) const {
    auto diameter = 1.0 / (1<<level);
    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0u; d<D; ++d)
        result[d]= { m_coords[cell][d]*diameter, (m_coords[cell][d]+1)*diameter };
    return result;
}

template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellForPoint(const std::array<double, D>& position, unsigned int targetLevel) {
    auto diameter = static_cast<double>(1 << targetLevel);

    std::array<uint32_t, D> coords;
    for (auto d = 0u; d < D; ++d)
        coords[d] = static_cast<uint32_t>(position[d] * diameter);

    return BitInterleaving<D>::interleave(coords, targetLevel);
}

template<unsigned int D>
bool SpatialTreeCoordinateHelper<D>::touching(unsigned int cellA, unsigned int cellB, unsigned int level) const {
    auto& coordA = m_coords[cellA];
    auto& coordB = m_coords[cellB];
    auto touching = true;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(coordA[d] - coordB[d]);
        dist = std::min(dist, (1<<level) - dist);
        touching &= dist<=1;
    }
    return touching;
}

template<unsigned int D>
double SpatialTreeCoordinateHelper<D>::dist(std::vector<double> &a, std::vector<double> &b) {
    assert(a.size() == b.size());
    assert(a.size() == D);

    // max over the torus distance in all dimensions
    auto result = 0.0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(a[d] - b[d]);
        dist = std::min(dist, 1.0-dist);
        result = std::max(result, dist);
    }
    return result;
}

template<unsigned int D>
double SpatialTreeCoordinateHelper<D>::dist(unsigned int cellA, unsigned int cellB, unsigned int level) const {

    // first work with integer d dimensional index
    auto& coordA = m_coords[cellA];
    auto& coordB = m_coords[cellB];
    auto result = 0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(coordA[d] - coordB[d]);
        dist = std::min(dist, (1<<level) - dist);
        result = std::max(result, dist);
    }

    // then apply the diameter
    auto diameter = 1.0 / (1<<level);
    return std::max(0.0, (result-1) * diameter); // TODO if cellA and cellB are not touching, this max is irrelevant
}


} // namespace girgs
