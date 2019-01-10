
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
unsigned int SpatialTreeCoordinateHelper<D>::cellForPoint(std::vector<double>& point, unsigned int targetLevel) const {
    // calculate coords
    assert(point.size() == D);
    auto diameter = static_cast<double>(1 << targetLevel);

    std::array<int, D> coords;
	for (auto d = 0u; d < D; ++d)
        coords[d] = static_cast<unsigned int>(point[d] * diameter);

	/*
	 * We now interleave the bits of the coordinates, let X[i,j] be the i-th bit (counting from LSB) of coordinate j,
	 * and let D' = D-1, and t = targetLevel-1. Then return
	 *  X[t, D'] o X[t, D'-1] o ... o X[t, 0]   o  X[t-1, D'] o ... o X[t-1, 0]   o  ...  o  X[0, D'] o ... o X[0, 0]
	 *
	 * The following in a naive implementation of this bit interleaving:
	 */
	unsigned int result = 0u;
	unsigned int bit = 0;
	for(auto l = 0u; l != targetLevel; l++) {
	    for(auto d = 0u; d != D; d++) {
	        result |= ((coords[d] >> l) & 1) << bit++;
	    }
	}

    return result + firstCellOfLevel(targetLevel);

	// TODO: We can use shifts and masks to use word-parallelism
	// TODO: We can use the AVX2 instruction pdep
	// TODO: Write benchmark what's faster given that we consider at most targetLevel bits per coordinate
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
