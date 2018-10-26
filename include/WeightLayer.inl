


template<unsigned int D>
WeightLayer<D>::WeightLayer(unsigned int layer,
                            unsigned int logW,
                            const SpatialTreeCoordinateHelper<D>& helper,
                            std::vector<Node<D>*>&& nodes)
    : m_layer(layer)
    , m_target_level((logW - (layer + 1)) / D) // w0*wi/W = 2^(-dl) solved for l
    , m_nodes(std::move(nodes))
{

    // a lot of assertions that we have the correct insertion level
    {
        // use layer+1 because paper has one based indices for weight layers
        auto v = std::pow(2,layer+1)/(1<<logW);
        auto volume_current = std::pow(2.0, -(m_target_level+0.0)*dimension);
        auto volume_next = std::pow(2.0, -(m_target_level+1.0)*dimension);
        auto deepestLevel = (logW-1)/dimension;
        assert(layer<logW); // because 2^ld = 1/W
        assert(v <= volume_current);
        assert(v >  volume_next);
        assert(m_target_level <= deepestLevel); // obvious
        assert(m_target_level >= 0);
    }

    // convenience constants
    const auto firstCell = firstCellOfLevel(m_target_level);
    const auto cellsInLevel = numCellsInLevel(m_target_level);
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
    for(auto i = 1; i < cellsInLevel; ++i){
        m_prefix_sums[i] = m_prefix_sums[i-1] + m_points_in_cell[i-1];
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
    assert(level <= m_target_level);
    assert(firstCellOfLevel(level) <= cell && cell < firstCellOfLevel(level+1)); // cell is from fromLevel

    // we want the begin-th and end-th cell in level targetLevel to be the first and last descendant of cell in this level
    // we could apply the firstChild function to find the first descendant but this is in O(1)
    auto descendants = numCellsInLevel(m_target_level - level);
    auto localIndexCell = cell - firstCellOfLevel(level);
    auto localIndexDescendant = localIndexCell * descendants; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendent
    auto begin = localIndexDescendant;
    auto end = begin + descendants - 1;

    assert(begin + firstCellOfLevel(level) < firstCellOfLevel(m_target_level+1));
    assert(end + firstCellOfLevel(level) < firstCellOfLevel(m_target_level+1));

    return m_prefix_sums[end] - m_prefix_sums[begin] + m_points_in_cell[end];
}


template<unsigned int D>
Node<D>* WeightLayer<D>::kthPoint(unsigned int cell, unsigned int level, int k) const {
    assert(level <= m_target_level);
    assert(firstCellOfLevel(level) <= cell && cell < firstCellOfLevel(level+1)); // cell is from fromLevel

    // same as in "pointsInCell"
    auto descendants = numCellsInLevel(m_target_level - level);
    auto localIndexCell = cell - firstCellOfLevel(level);
    auto localIndexDescendant = localIndexCell * descendants;
    auto begin = localIndexDescendant;

    return m_A[m_prefix_sums[begin] + k];
}

