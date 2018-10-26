


template<unsigned int D>
std::vector<Node<D>> SpatialTree<D>::generateGraph(std::vector<double>& weights) {

    auto graph = Graph(weights.size());

    // sample positions
    std::mt19937 gen; // TODO expose seed
    std::uniform_real_distribution<> dist;
    for(auto i=0; i<weights.size(); ++i) {
        auto &node = graph[i];
        node.weight = weights[i];
        node.index = i;
        for (auto d = 0u; d < D; ++d)
            node.coord[d] = dist(gen);
    }

    // determine size of tree and weight layers
    m_W = std::accumulate(weights.begin(), weights.end(), 0.0);
    m_layers = static_cast<unsigned int>(std::log2(m_W));
    m_levels = (m_layers - 1) / D + 1; // level are in range [0, (layer-1)/D]
    m_helper = SpatialTreeCoordinateHelper<D>(m_levels);

    // sort weights into exponentially growing layers
    {   // block to let weightLayerNodes go out of scope after it was moved away
        auto weightLayerNodes = std::vector<std::vector<Node<D> *>>(m_layers);
        for (auto i = 0; i < weights.size(); ++i)
            weightLayerNodes[std::log2(weights[i])].push_back(&graph[i]);

        // build spatial structure described in paper
        for (auto layer = 0; layer < m_layers; ++layer)
            m_weight_layers.emplace_back(layer, m_layers, m_helper, std::move(weightLayerNodes[layer]));
    }

    // sample all edges
    visitCellPair(0,0,0, graph);

    return graph;
}


template<unsigned int D>
void SpatialTree<D>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, std::vector<Node<D>> &graph) const {

    auto touching = m_helper.touching(cellA, cellB, level);

    if(cellA == cellB or touching) {

        // sample all type 1 occurrences with this cell pair
        // TODO dont do bruteforce over all layer-pairs
        for(auto i=0; i<m_layers; ++i)
            for (auto j=0; j<m_layers; ++j)
                if(std::max(((int)m_layers-(i+j+2))/(int)dimension, 0) == level)
                    sampleBetweenViAVjB(cellA, cellB, level, i, j, graph);

    } else { // not touching

        // sample all type 2 occurrences with this cell pair
        for(auto i=0; i<m_layers; ++i)
            for (auto j=0; j<m_layers; ++j)
                if(std::max(((int)m_layers-(i+j+2))/(int)dimension, 0) >= level){
                    sampleBetweenViAVjB(cellA, cellB, level, i, j, graph);
                } else {
                    break; // if condition failed it will also fail for all subsequent j
                }
    }

    // break if last level reached
    if(level >= (m_layers-2)/dimension) // we skip level (L - 1)/d, because i,j are compared at most in depth (L-2)/d
        return;

    if(cellA == cellB){
        // recursive call for all children pairs of cell A=B
        for(auto a = firstChild(cellA); a<firstChild(cellA)+numChildren; ++a)
            for(auto b = a; b<firstChild(cellA)+numChildren; ++b)
                visitCellPair(a,b, level+1, graph);
    }


    if(touching && cellA!=cellB) {
        // recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they don't
        for(auto a = firstChild(cellA); a<firstChild(cellA)+numChildren; ++a)
            for(auto b = firstChild(cellB); b<firstChild(cellB)+numChildren; ++b)
                visitCellPair(a, b, level+1, graph);
    }
}


template<unsigned int D>
void SpatialTree<D>::sampleBetweenViAVjB(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j, SpatialTree::Graph &graph) const
{
    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =0; kB<sizeV_j_B; ++kB) {
            auto nodeInA = m_weight_layers[i].kthPoint(cellA, level, kA);
            auto nodeInB = m_weight_layers[j].kthPoint(cellB, level, kB);

            // points are in correct cells
            assert(cellA == m_helper.cellForPoint(nodeInA->coord, level));
            assert(cellB == m_helper.cellForPoint(nodeInB->coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA->weight)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB->weight)));

            // TODO replace by real sampling
            nodeInA->edges.push_back(nodeInB);
            if(cellA != cellB) // to sample once even if ViA = ViB
                nodeInB->edges.push_back(nodeInA);
        }
    }
}
