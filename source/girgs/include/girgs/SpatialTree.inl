

namespace girgs {


template<unsigned int D>
void SpatialTree<D>::generateEdges(std::vector<Node>& graph, double alpha, int seed) {

    // init member and determine min max and sum of weights
    m_alpha = alpha;
    m_w0 = std::numeric_limits<double>::infinity();
    m_wn = 0.0;
    m_W = 0.0;
    for(auto i=0u; i<graph.size(); ++i) {
        graph[i].index = i;
        m_w0 = std::min(m_w0, graph[i].weight);
        m_wn = std::max(m_wn, graph[i].weight);
        m_W += graph[i].weight;
    }

    // determine size of tree and weight layers
    m_layers = static_cast<unsigned int>(floor(std::log2(m_wn/m_w0)))+1;
    m_baseLevelConstant = static_cast<int>(std::log2(m_W/m_w0/m_w0)); // log2(W/w0^2)
    m_levels = partitioningBaseLevel(0,0) + 1; // (log2(W/w0^2) - 2) / d
    // we need a helper for the deepest insertion level which possibly is one larger than the deepest comparison level
    m_helper = SpatialTreeCoordinateHelper<D>(m_levels+1);

    // determine which layer pairs to sample in which level
    // TODO maybe also save type2 pairs? or loop over multiple vectors for type2?
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(i, j)].emplace_back(i,j);


    // sort weights into exponentially growing layers
    {   // block to let weightLayerNodes go out of scope after it was moved away
        auto weightLayerNodes = std::vector<std::vector<Node*>>(m_layers);
        for (auto i = 0u; i < graph.size(); ++i)
            weightLayerNodes[std::log2(graph[i].weight/m_w0)].push_back(&graph[i]);

        // build spatial structure described in paper
        for (auto layer = 0u; layer < m_layers; ++layer)
            m_weight_layers.emplace_back(layer, weightLayerTargetLevel(layer), m_helper, std::move(weightLayerNodes[layer]));
    }

    // one random generator and distribution for each thread
    const auto num_threads = omp_get_max_threads();
    m_gens.resize(num_threads);
    m_dists.resize(num_threads);
    for (int thread = 0; thread < num_threads; thread++) {
        m_gens[thread].seed(seed >= 0 ? seed+thread : std::random_device()());
    } 

#ifndef NDEBUG
    // ensure that all node pairs are compared either type 1 or type 2
    m_type1_checks.resize(num_threads, 0);
    m_type2_checks.resize(num_threads, 0);
#endif // NDEBUG

    // sample all edges
	if (num_threads == 1) { 
        // sequential
		visitCellPair(0, 0, 0);
    } else {
        // parallel see docs for visitCellPair_sequentialStart
        const auto first_parallel_level = static_cast<unsigned int>(std::ceil(std::log2(4.0*num_threads) / D));
        const auto parallel_cells = SpatialTreeCoordinateHelper<D>::numCellsInLevel(first_parallel_level);
        const auto first_parallel_cell = SpatialTreeCoordinateHelper<D>::firstCellOfLevel(first_parallel_level);

        // saw off recursion before "first_parallel_level" and save all calls that would be made 
        auto parallel_calls = std::vector<std::vector<unsigned int>>(parallel_cells);
        visitCellPair_sequentialStart(0, 0, 0, first_parallel_level, parallel_calls);
        
        // do the collected calls in parallel
        #pragma omp parallel for schedule(static), num_threads(num_threads) // dynamic scheduling would be better but not reproducible
        for (int i = 0; i < parallel_cells; ++i) {
            auto current_cell = first_parallel_cell + i;
            for (auto each : parallel_calls[i])
                visitCellPair(current_cell, each, first_parallel_level);
        }
    }

#ifndef NDEBUG
    // after sampling the graph the sum of type 1 and type 2 checks should always be n(n-1) to ensure that all edges were considered
    auto type1 = std::accumulate(m_type1_checks.begin(), m_type1_checks.end(), 0ll);
    auto type2 = std::accumulate(m_type2_checks.begin(), m_type2_checks.end(), 0ll);
    assert(type1 + type2 == graph.size()*(graph.size() - 1ll)
        || alpha == std::numeric_limits<double>::infinity()); // we do not compare all nodes in threshold since we skip all type 2 checks
#endif // NDEBUG 
}


template<unsigned int D>
void SpatialTree<D>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {
    using Helper = SpatialTreeCoordinateHelper<D>;

    auto touching = m_helper.touching(cellA, cellB, level);
    if(cellA == cellB || touching) {

        // sample all type 1 occurrences with this cell pair
        for(auto& layer_pair : m_layer_pairs[level]){
            assert(partitioningBaseLevel(layer_pair.first, layer_pair.second) == level);
            if(cellA != cellB || layer_pair.first <= layer_pair.second)
                sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
        }

    } else { // not touching
		if (m_alpha == std::numeric_limits<double>::infinity())
			return;
        // sample all type 2 occurrences with this cell pair
        // TODO is it correct to sample all pairs in m_layer_pairs[k] for k in range(level, m_levels) ???
        for(auto i=0u; i<m_layers; ++i)
            for (auto j=0u; j<m_layers; ++j)
                if(partitioningBaseLevel(i,j) >= level){
                    sampleTypeII(cellA, cellB, level, i, j);
                } else {
                    break; // if condition failed it will also fail for all subsequent j
                }
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    if(touching) {
        // recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they don't
        for(auto a = Helper::firstChild(cellA); a<=Helper::lastChild(cellA); ++a)
            for(auto b = cellA == cellB ? a : Helper::firstChild(cellB); b<=Helper::lastChild(cellB); ++b)
                visitCellPair(a, b, level+1);
    }
}



template<unsigned int D>
void SpatialTree<D>::visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
                                                   unsigned int first_parallel_level,
                                                   std::vector<std::vector<unsigned int>> &parallel_calls) {
    using Helper = SpatialTreeCoordinateHelper<D>;

    auto touching = m_helper.touching(cellA, cellB, level);
    if(cellA == cellB || touching) {
        // sample all type 1 occurrences with this cell pair
        for(auto& layer_pair : m_layer_pairs[level]){
            assert(partitioningBaseLevel(layer_pair.first, layer_pair.second) == level);
            if(cellA != cellB || layer_pair.first <= layer_pair.second)
                sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
        }
    } else { // not touching
        // sample all type 2 occurrences with this cell pair
        for(auto i=0u; i<m_layers; ++i)
            for (auto j=0u; j<m_layers; ++j)
                if(partitioningBaseLevel(i,j) >= level){
                    sampleTypeII(cellA, cellB, level, i, j);
                } else {
                    break; // if condition failed it will also fail for all subsequent j
                }
    }

    // break if last level reached
    if (level == m_levels - 1) // if we are at the last level we don't need recursive calls
        return;

    if(touching) {
        // recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they don't
        for(auto a = Helper::firstChild(cellA); a<=Helper::lastChild(cellA); ++a)
            for(auto b = cellA == cellB ? a : Helper::firstChild(cellB); b<=Helper::lastChild(cellB); ++b){
                if(level+1 == first_parallel_level)
                    parallel_calls[a-Helper::firstCellOfLevel(first_parallel_level)].push_back(b);
                else
                    visitCellPair_sequentialStart(a, b, level+1, first_parallel_level, parallel_calls);
            }

    }
}



template<unsigned int D>
void SpatialTree<D>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{

    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
	if (sizeV_i_A == 0 || sizeV_j_B == 0)
		return;

#ifndef NDEBUG
    m_type1_checks[omp_get_thread_num()] += (cellA == cellB && i == j) 
        ? sizeV_i_A * (sizeV_i_A-1)  // all pairs in AxA without {v,v}
        : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
#endif // NDEBUG

	Node* const * firstA = m_weight_layers[i].firstPointPointer(cellA, level);
	Node* const * firstB = m_weight_layers[j].firstPointPointer(cellB, level);

    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =(cellA == cellB && i==j ? kA+1 : 0); kB<sizeV_j_B; ++kB) {
            Node* nodeInA = *(firstA+kA);
            Node* nodeInB = *(firstB+kB);
			//Node* nodeInA = m_weight_layers[i].kthPoint(cellA, level, kA);
			//Node* nodeInB = m_weight_layers[j].kthPoint(cellB, level, kB);

			// pointer magic gives same results
			assert(nodeInA == m_weight_layers[i].kthPoint(cellA, level, kA));
			assert(nodeInB == m_weight_layers[j].kthPoint(cellB, level, kB));

            // points are in correct cells
            assert(cellA == m_helper.cellForPoint(nodeInA->coord, level));
            assert(cellB == m_helper.cellForPoint(nodeInB->coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA->weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB->weight/m_w0)));

            assert(nodeInA->index != nodeInB->index);
            auto dist = m_helper.dist(nodeInA->coord, nodeInB->coord);
            if(checkEdgeExplicit(dist, nodeInA->weight, nodeInB->weight)){
                nodeInA->edges.push_back(nodeInB);
                //nodeInB->edges.push_back(nodeInA);
            }
        }
    }
}


template<unsigned int D>
void SpatialTree<D>::sampleTypeII(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    long long sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    long long sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
    if(sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

    // get upper bound for probability
    auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    auto dist_lower_bound = std::pow(m_helper.dist(cellA, cellB, level), dimension);
    auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    if(max_connection_prob <= 1e-10)
        return;

    // if we must sample all pairs we treat this as type 1 sampling
    // also, 1.0 is no valid prob for a geometric dist (see c++ std)
    if(max_connection_prob == 1.0){
        sampleTypeI(cellA, cellB, level, i, j);
        return;
    }

    // init geometric distribution
    auto threadID = omp_get_thread_num();
    auto& gen = m_gens[threadID];
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);

#ifndef NDEBUG
    m_type2_checks[omp_get_thread_num()] += 2 * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    for (auto r = geo(gen); r < sizeV_i_A * sizeV_j_B; r += 1 + geo(gen)) {
        // determine the r-th pair
        Node* nodeInA = m_weight_layers[i].kthPoint(cellA, level, r%sizeV_i_A);
        Node* nodeInB = m_weight_layers[j].kthPoint(cellB, level, r/sizeV_i_A);

        // points are in correct cells
        assert(cellA == m_helper.cellForPoint(nodeInA->coord, level));
        assert(cellB == m_helper.cellForPoint(nodeInB->coord, level));

        // points are in correct weight layer
        assert(i == static_cast<unsigned int>(std::log2(nodeInA->weight/m_w0)));
        assert(j == static_cast<unsigned int>(std::log2(nodeInB->weight/m_w0)));

        // get actual connection probability
        auto w = nodeInA->weight*nodeInB->weight/m_W;
        auto d = std::pow(m_helper.dist(nodeInA->coord, nodeInB->coord), dimension);
        auto connection_prob = std::min(std::pow(w/d, m_alpha), 1.0);
        assert(w < w_upper_bound);
        assert(d >= dist_lower_bound);

        if(m_dists[threadID](gen) < connection_prob/max_connection_prob) {
            nodeInA->edges.push_back(nodeInB);
            //nodeInB->edges.push_back(nodeInA);
        }
    }
}


template<unsigned int D>
unsigned int SpatialTree<D>::weightLayerTargetLevel(int layer) const {
    // -1 coz w0 is the upper bound for layer 0 in paper and our layers are shifted by -1
    auto result = std::max((m_baseLevelConstant - layer - 1) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct insertion level
        assert(0 <= layer && layer < m_layers);
        assert(0 <= result && result <= m_levels); // note the result may be one larger than the deepest level (hence the <= m_levels)
        auto volume_requested  = m_w0*m_w0*std::pow(2,layer+1)/m_W; // v(i) = w0*wi/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


template<unsigned int D>
unsigned int SpatialTree<D>::partitioningBaseLevel(int layer1, int layer2) const {

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct comparison level
        assert(0 <= layer1 && layer1 < m_layers);
        assert(0 <= layer2 && layer2 < m_layers);
        assert(0 <= result && result < m_levels || m_levels == 0);
        auto volume_requested  = m_w0*std::pow(2,layer1+1) * m_w0*std::pow(2,layer2+1) / m_W; // v(i,j) = wi*wj/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


template<unsigned int D>
bool SpatialTree<D>::checkEdgeExplicit(double dist, double w1, double w2) {
    auto w_term = w1*w2/m_W;
	auto d_term = 1.0; // dist^D
	for (int i = 0; i < D; ++i)
		d_term *= dist;
    if(m_alpha == std::numeric_limits<double>::infinity())
        return d_term < w_term;

    auto edge_prob = std::min(std::pow(w_term/d_term, m_alpha), 1.0);
    auto threadID = omp_get_thread_num();
    return m_dists[threadID](m_gens[threadID]) < edge_prob;
}



} // namespace girgs
