


template<unsigned int D>
std::vector<Node> SpatialTree<D>::generateGraph(const std::vector<double> &weights, double alpha, double c, int seed) {

    // init member
    m_alpha = alpha;
    m_c = c;
    m_1 = 0;
    m_2 = 0;
    m_gen = std::mt19937(seed >= 0 ?  seed : std::random_device()());
    m_dist = std::uniform_real_distribution<>(0.0, 1.0);

    // determine size of tree and weight layers
    m_W = std::accumulate(weights.begin(), weights.end(), 0.0);
    m_layers = static_cast<unsigned int>(std::log2(m_W));
    m_levels = (m_layers - 1) / D + 1; // level are in range [0, (layer-1)/D]
    m_helper = SpatialTreeCoordinateHelper<D>(m_levels);

    // sample positions
    auto graph = Graph(weights.size());
    for(auto i=0; i<weights.size(); ++i) {
        auto &node = graph[i];
        node.weight = weights[i];
        node.index = i;
        for (auto d = 0u; d < D; ++d)
            node.coord.push_back(m_dist(m_gen));
    }

    // sort weights into exponentially growing layers
    {   // block to let weightLayerNodes go out of scope after it was moved away
        auto weightLayerNodes = std::vector<std::vector<Node*>>(m_layers);
        for (auto i = 0; i < weights.size(); ++i)
            weightLayerNodes[std::log2(weights[i])].push_back(&graph[i]);

        // build spatial structure described in paper
        for (auto layer = 0; layer < m_layers; ++layer)
            m_weight_layers.emplace_back(layer, m_layers, m_helper, std::move(weightLayerNodes[layer]));
    }

    // sample all edges
    visitCellPair(0,0,0, graph);

    std::cout << m_1 << "\t" << m_2 << std::endl;
    assert(m_1 + m_2 == graph.size() * graph.size());
    return graph;
}


template<unsigned int D>
void SpatialTree<D>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, std::vector<Node> &graph) {

    auto touching = m_helper.touching(cellA, cellB, level);

    // TODO early return if A or B empty

    if(cellA == cellB or touching) {

        // sample all type 1 occurrences with this cell pair
        // TODO dont do bruteforce over all layer-pairs
        for(auto i=0; i<m_layers; ++i)
            for (auto j=(cellA == cellB ? i : 0); j<m_layers; ++j)
                if(std::max(((int)m_layers-(i+j+2))/(int)dimension, 0) == level)
                    sampleTypeI(cellA, cellB, level, i, j, graph);

    } else { // not touching

        // sample all type 2 occurrences with this cell pair
        for(auto i=0; i<m_layers; ++i)
            for (auto j=0; j<m_layers; ++j)
                if(std::max(((int)m_layers-(i+j+2))/(int)dimension, 0) >= level){
                    sampleTypeII(cellA, cellB, level, i, j, graph);
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
void SpatialTree<D>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j, SpatialTree::Graph &graph)
{

    // a lot of assertions that we have the correct comparison level
    {
        // use layer+1 because paper has one based indices for weight layers
        auto v = std::pow(2.0,i+1.0) * std::pow(2.0,j+1.0) / (m_W);             // in paper v(i,j)
        auto volume_comparison_level = std::pow(2.0, -(level+0.0)*dimension);   // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper       = std::pow(2.0, -(level+1.0)*dimension);
        assert(level < m_levels);
        assert(level >= 0);
        assert(v <= volume_comparison_level || v >= 1.0);
        assert(v >  volume_one_deeper);
    }

    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);

    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =(cellA == cellB && i==j ? kA : 0); kB<sizeV_j_B; ++kB) {
            Node* nodeInA = m_weight_layers[i].kthPoint(cellA, level, kA);
            Node* nodeInB = m_weight_layers[j].kthPoint(cellB, level, kB);

            // points are in correct cells
            assert(cellA == m_helper.cellForPoint(nodeInA->coord, level));
            assert(cellB == m_helper.cellForPoint(nodeInB->coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA->weight)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB->weight)));

            // TODO replace by real sampling
            m_1 += 1 + (nodeInA->index != nodeInB->index);
            auto dist = m_helper.dist(nodeInA->coord, nodeInB->coord);
            if(checkEdgeExplicit(dist, nodeInA->weight, nodeInB->weight)){
                nodeInA->edges.push_back(nodeInB);
                if(nodeInA->index != nodeInB->index)
                    nodeInB->edges.push_back(nodeInA);
            }
        }
    }
}


template<unsigned int D>
void SpatialTree<D>::sampleTypeII(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j, SpatialTree::Graph &graph)
{

    // a lot of assertions that we have the correct comparison level
    {
        // use layer+1 because paper has one based indices for weight layers
        auto v = std::pow(2.0,i+1.0) * std::pow(2.0,j+1.0) / (m_W);             // in paper v(i,j)
        auto volume_comparison_level = std::pow(2.0, -(level+0.0)*dimension);   // in paper \mu with v <= \mu < O(v)
        assert(level < m_levels);
        assert(level >= 0);
        assert(v <= volume_comparison_level || v >= 1.0);
        assert(cellA != cellB);
        // the level of type 1 pairs in the partitioning we belong to
        auto partitioning_base_level = ((int)m_layers-(i+j+2))/(int)dimension;
        assert(0 < partitioning_base_level);
        assert(level <= partitioning_base_level);
    }

    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
    if(sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

    m_2 += 2*sizeV_i_A*sizeV_j_B;

    // implicit sampling
    auto w_upper_bound = (1<<i+1) * (1<<j+1) / m_W;
    auto dist_lower_bound = std::pow(m_helper.dist(cellA, cellB, level), dimension);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    auto max_connection_prob = std::min(m_c*std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    auto geo = [this](double p) -> int {
        auto R = this->m_dist(this->m_gen);
        return p==1 ? 1 : std::ceil(std::log2(R) / std::log2(1-p)); // this does not work if p=1
    };
    auto r = geo(max_connection_prob);
    while(r <= sizeV_i_A * sizeV_j_B){
        // determine the r-th pair
        Node *nodeInA = m_weight_layers[i].kthPoint(cellA, level, (r-1)%sizeV_i_A);
        Node *nodeInB = m_weight_layers[j].kthPoint(cellB, level, (r-1)/sizeV_i_A);

        auto w = nodeInA->weight*nodeInB->weight/m_W;
        assert(w < w_upper_bound);
        auto d = std::pow(m_helper.dist(nodeInA->coord, nodeInB->coord), dimension);
        assert(d >= dist_lower_bound);
        auto connection_prob = std::min(m_c*std::pow(w/d, m_alpha), 1.0);
        if(m_dist(m_gen) < connection_prob/max_connection_prob) {
            nodeInA->edges.push_back(nodeInB);
            nodeInB->edges.push_back(nodeInA);
        }
        r += geo(max_connection_prob);
    }
}



template<unsigned int D>
bool SpatialTree<D>::checkEdgeExplicit(double dist, double w1, double w2) {

    double p;
    // TODO do this check earlier
    if(m_alpha == std::numeric_limits<double>::infinity()){
        p = dist < m_c*std::pow(w1*w2/m_W, 1.0/D);
    } else {
        auto w = std::pow(w1*w2/m_W, m_alpha);
        auto d = std::pow(dist, m_alpha * dimension);
        p = std::min(m_c*w/d, 1.0);
    }

    return m_dist(m_gen) < p;
}

