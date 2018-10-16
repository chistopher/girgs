


template<unsigned int D, unsigned int L>
constexpr SpatialTree<D,L>::SpatialTree()
    : m_coords()
    , m_coords2Index()
    , m_points_in_cell()
    , m_prefix_sums()
    , m_A(nullptr)
{

    // calculate coords // TODO find implicit way to get from index to coords and back
    for(auto d=0; d<D; ++d) m_coords[0][d] = 0; // first cell has [0]^D // no "m_coords[0].fill(0);" because that is not constexpr

    for(auto l=1; l<L; ++l) {
        for(auto cell = firstCellOfLevel(l); cell < firstCellOfLevel(l+1); ++cell) {
            m_coords[cell] = m_coords[parent(cell)];
            auto bitmask = ((cell-1)&(numChildren-1)); // & is unnecessary but removes all but the last D bits
            for(auto d=0; d<D; ++d) { // check orientation to parent in all dimensions with d-th bit in bitmask
                m_coords[cell][d] *= 2; // every level doubles each dimension
                m_coords[cell][d] += static_cast<bool>(bitmask&(1<<d)); // add 1 if we are the "right" of parent in dimension d
            }

            // save inverse mapping
            auto packedCoord = 0;
            for(auto d=0; d<D; ++d)
                packedCoord += m_coords[cell][d] * (1<<(l*d));
            m_coords2Index[firstCellOfLevel(l) + packedCoord] = cell;
        }
    }
}


template<unsigned int D, unsigned int L>
int SpatialTree<D,L>::pointsInCell(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel) const {
    assert(fromLevel <= targetLevel);
    assert(firstCellOfLevel(fromLevel) <= cell && cell < firstCellOfLevel(fromLevel+1)); // cell is from fromLevel

    // we want the begin-th and end-th cell in level targetLevel to be the first and last descendent of cell in this level
    // we could apply the firstChild function to find the first descendent but this is in O(1)
    auto descendents = numCellsInLevel(targetLevel - fromLevel); // 2^(D*L)
    auto localIndexCell = cell - firstCellOfLevel(fromLevel);
    auto localIndexDescendent = localIndexCell * descendents; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendent
    auto begin = localIndexDescendent + firstCellOfLevel(targetLevel);
    auto end = begin + descendents - 1;

    assert(begin < firstCellOfLevel(targetLevel+1));
    assert(end < firstCellOfLevel(targetLevel+1));

    return m_prefix_sums[end] - m_prefix_sums[begin] + m_points_in_cell[end];
}


template<unsigned int D, unsigned int L>
std::array<std::pair<double, double>, D> SpatialTree<D, L>::bounds(int cell, unsigned int level) const {
    auto diameter = 1.0 / (1<<level);
    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0; d<D; ++d)
        result[d]= { m_coords[cell][d]*diameter, (m_coords[cell][d]+1)*diameter };
    return result;
}

template<unsigned int D, unsigned int L>
unsigned int SpatialTree<D, L>::cellForPoint(std::array<double, D>& point, unsigned int targetLevel) const {

    // calculate coords
    auto diameter = 1.0 / (1<<targetLevel);
    std::array<int, D> coords;
    for(auto d=0; d<D; ++d) {
        coords[d] = static_cast<int>(point[d] / diameter);
        assert(coords[d] < (1<<targetLevel));
    }

    // get from coords to index
    auto packedCoord = 0;
    for(auto d=0; d<D; ++d)
        packedCoord += coords[d] * (1<<(targetLevel*d));
    return m_coords2Index[firstCellOfLevel(targetLevel) + packedCoord];
}
