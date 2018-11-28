
namespace girgs {

template<unsigned int D>
WeightLayer<D>::WeightLayer(unsigned int layer,
                            unsigned int targetLevel,
                            const SpatialTreeCoordinateHelper<D>& helper,
                            std::vector<Node*>&& nodes)
    : m_layer(layer)
    , m_target_level(targetLevel) // w0*wi/W = 2^(-dl) solved for l --- l = (log2(W/w0^2) - i) / d
    , m_nodes(std::move(nodes))
{

    // convenience constants
    const auto firstCell = SpatialTreeCoordinateHelper<D>::firstCellOfLevel(m_target_level);
    const auto cellsInLevel = SpatialTreeCoordinateHelper<D>::numCellsInLevel(m_target_level);
    const auto lastCell = firstCell + cellsInLevel - 1;

    // allocate stuff
    m_points_in_cell.resize(cellsInLevel, 0);
    m_prefix_sums.resize(cellsInLevel, 0);
    m_A.resize(m_nodes.size(), nullptr);

    // count num of points in each cell
    for(auto node : m_nodes){
        auto targetCell = helper.cellForPoint(node->coord, m_target_level);
        assert(firstCell <= targetCell && targetCell <= lastCell); // cell on right level
        ++m_points_in_cell[targetCell-firstCell];
    }

    // fill prefix sums
    // prefix_sums[i] is the number of all points in cells j<i of the same level
    m_prefix_sums[0] = 0;
    for(auto cell = 1u; cell < cellsInLevel; ++cell) {
        m_prefix_sums[cell] = m_prefix_sums[cell-1] + m_points_in_cell[cell-1];
    }

    // fill point lookup
    auto num_inserted = std::vector<int>(cellsInLevel, 0); // keeps track of bucket size for counting sort
    for(auto node : m_nodes){
        auto targetCell = helper.cellForPoint(node->coord, m_target_level);
        assert(firstCell <= targetCell && targetCell <= lastCell); // cell on right level
        m_A[m_prefix_sums[targetCell-firstCell] + num_inserted[targetCell-firstCell]] = node;
        ++num_inserted[targetCell-firstCell];
    }
}


template<unsigned int D>
int WeightLayer<D>::pointsInCell(unsigned int cell, unsigned int level) const {
    using Helper = SpatialTreeCoordinateHelper<D>;
    assert(level <= m_target_level);
    assert(Helper::firstCellOfLevel(level) <= cell && cell < Helper::firstCellOfLevel(level+1)); // cell is from correct level

    // we want the begin-th and end-th cell in level targetLevel to be the first and last descendant of cell in this level
    // we could apply the firstChild function to find the first descendant but this is in O(1)
    auto descendants = Helper::numCellsInLevel(m_target_level - level);
    auto localIndexCell = cell - Helper::firstCellOfLevel(level);
    auto localIndexDescendant = localIndexCell * descendants; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendant
    auto begin = localIndexDescendant;
    auto end = begin + descendants - 1;

    assert(begin + Helper::firstCellOfLevel(level) < Helper::firstCellOfLevel(m_target_level+1));
    assert(end + Helper::firstCellOfLevel(level) < Helper::firstCellOfLevel(m_target_level+1));

    return m_prefix_sums[end] - m_prefix_sums[begin] + m_points_in_cell[end];
}


template<unsigned int D>
Node* WeightLayer<D>::kthPoint(unsigned int cell, unsigned int level, int k) const {
    using Helper = SpatialTreeCoordinateHelper<D>;
    assert(level <= m_target_level);
    assert(Helper::firstCellOfLevel(level) <= cell && cell < Helper::firstCellOfLevel(level+1)); // cell is from fromLevel

    // same as in "pointsInCell"
    auto descendants = Helper::numCellsInLevel(m_target_level - level);
    auto localIndexCell = cell - Helper::firstCellOfLevel(level);
    auto localIndexDescendant = localIndexCell * descendants;
    auto begin = localIndexDescendant;

    return m_A[m_prefix_sums[begin] + k];
}


} // namespace girgs
