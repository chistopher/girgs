#include <girgs/ScopedTimer.h>
#include <girgs/Helper.h>

namespace girgs {


template<unsigned int D, typename EdgeCallback>
SpatialTree<D, EdgeCallback>::SpatialTree(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, double alpha, EdgeCallback& edgeCallback, bool profile)
: m_EdgeCallback(edgeCallback)
, m_profile(profile)
, m_alpha(alpha)
, m_n(weights.size())
, m_w0(*std::min_element(weights.begin(), weights.end()))
, m_wn(*std::max_element(weights.begin(), weights.end()))
, m_W(std::accumulate(weights.begin(), weights.end(), 0.0))
, m_baseLevelConstant(static_cast<int>(std::log2(m_W/m_w0/m_w0))) // log2(W/w0^2)
, m_layers(static_cast<unsigned int>(floor(std::log2(m_wn/m_w0)))+1)
, m_levels(partitioningBaseLevel(0,0) + 1) // (log2(W/w0^2) - 2) / d
, m_helper(SpatialTreeCoordinateHelper<D>(m_levels+1)) // helper for the deepest insertion level which possibly is one larger than the deepest comparison level
{
    assert(weights.size() == positions.size());
    assert(positions.size() > 0 && positions.front().size() == D);

    ScopedTimer timer("Preprocessing", profile);

    // determine which layer pairs to sample in which level
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(i, j)].emplace_back(i,j);

    // sort weights into exponentially growing layers
    {   // block to let weightLayerNodes go out of scope after it was moved away
        auto weightLayerNodes = std::vector<std::vector<int>>(m_layers);
        for (auto i = 0u; i < weights.size(); ++i)
            weightLayerNodes[std::log2(weights[i]/m_w0)].push_back(i);

        // build spatial structure described in paper
        for (auto layer = 0u; layer < m_layers; ++layer)
            m_weight_layers.push_back(
                    WeightLayer<D>(layer, weightLayerTargetLevel(layer), m_helper, std::move(weightLayerNodes[layer]), weights, positions)
            );
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::generateEdges(int seed) {

    // one random generator and distribution for each thread
    const auto num_threads = omp_get_max_threads();
    m_gens.resize(num_threads);
    for (int thread = 0; thread < num_threads; thread++) {
        m_gens[thread].seed(seed >= 0 ? seed+thread : std::random_device()());
    } 

#ifndef NDEBUG
    // ensure that all node pairs are compared either type 1 or type 2
    m_type1_checks = 0;
    m_type2_checks = 0;
#endif // NDEBUG

    // sample all edges
	if (num_threads == 1) { 
        // sequential
		visitCellPair(0, 0, 0);
		assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
		return;
    }

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

    assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {
    using Helper = SpatialTreeCoordinateHelper<D>;

    if(!m_helper.touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
		if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = Helper::firstChild(cellA); a<=Helper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : Helper::firstChild(cellB); b<=Helper::lastChild(cellB); ++b)
            visitCellPair(a, b, level+1);
}



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
                                                   unsigned int first_parallel_level,
                                                   std::vector<std::vector<unsigned int>> &parallel_calls) {
    using Helper = SpatialTreeCoordinateHelper<D>;

    if(!m_helper.touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
		if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if (level == m_levels - 1) // if we are at the last level we don't need recursive calls
        return;

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



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) == level
        || !m_helper.touching(cellA, cellB, level)); // in this case we were redirected from typeII with maxProb==1.0
    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
	if (sizeV_i_A == 0 || sizeV_j_B == 0)
		return;

#ifndef NDEBUG
    #pragma omp atomic
    m_type1_checks += (cellA == cellB && i == j)
            ? sizeV_i_A * (sizeV_i_A-1)  // all pairs in AxA without {v,v}
            : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
#endif // NDEBUG

    auto threadId = omp_get_thread_num();
    std::uniform_real_distribution<> dist;
    const auto* firstA = m_weight_layers[i].firstPointPointer(cellA, level);
    const auto* firstB = m_weight_layers[j].firstPointPointer(cellB, level);
    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =(cellA == cellB && i==j ? kA+1 : 0); kB<sizeV_j_B; ++kB) {
            const Node<D>& nodeInA = firstA[kA];
            const Node<D>& nodeInB = firstB[kB];

			// pointer magic gives same results
			assert(nodeInA.index == m_weight_layers[i].kthPoint(cellA, level, kA).index);
			assert(nodeInB.index == m_weight_layers[j].kthPoint(cellB, level, kB).index);

            // points are in correct cells
            assert(cellA - m_helper.firstCellOfLevel(level) == m_helper.cellForPoint({nodeInA.coord.begin(), nodeInA.coord.end()}, level));
            assert(cellB - m_helper.firstCellOfLevel(level) == m_helper.cellForPoint({nodeInB.coord.begin(), nodeInB.coord.end()}, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

            assert(nodeInA.index != nodeInB.index);
            auto distance = nodeInA.distance(nodeInB);
            auto w_term = nodeInA.weight*nodeInB.weight/m_W;
            const auto d_term = pow_to_the<D>(distance);

            if(m_alpha == std::numeric_limits<double>::infinity()) {
                if(d_term < w_term)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            } else {
                auto edge_prob = std::min(std::pow(w_term/d_term, m_alpha), 1.0);
                if(dist(m_gens[threadId]) < edge_prob)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            }
        }
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeII(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) >= level);
    long long sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    long long sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
    if(sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

    // get upper bound for probability
    const auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    const auto cell_distance = m_helper.dist(cellA, cellB, level);
    const auto dist_lower_bound = pow_to_the<D>(cell_distance);
    const auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    const auto num_pairs = sizeV_i_A * sizeV_j_B;
    const auto expected_samples = num_pairs * max_connection_prob;

    // if we must sample all pairs we treat this as type 1 sampling
    // also, 1.0 is no valid prob for a geometric dist (see c++ std)
    if(max_connection_prob == 1.0){
        sampleTypeI(cellA, cellB, level, i, j);
        return;
    }

#ifndef NDEBUG
    #pragma omp atomic
    m_type2_checks += 2 * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    if(expected_samples < 1e-6)
        return;

    // init geometric distribution
    auto threadId = omp_get_thread_num();
    auto& gen = m_gens[threadId];
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);
    std::uniform_real_distribution<> dist;
    const auto* firstA = m_weight_layers[i].firstPointPointer(cellA, level);
    const auto* firstB = m_weight_layers[j].firstPointPointer(cellB, level);
    for (auto r = geo(gen); r < num_pairs; r += 1 + geo(gen)) {
        // determine the r-th pair
        const Node<D>& nodeInA = firstA[r%sizeV_i_A];
        const Node<D>& nodeInB = firstB[r/sizeV_i_A];

        // points are in correct weight layer
        assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
        assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

        // get actual connection probability
        const auto distance = nodeInA.distance(nodeInB);
        const auto w_term = nodeInA.weight*nodeInB.weight/m_W;
        const auto d_term = pow_to_the<D>(distance);
        const auto connection_prob = std::min(std::pow(w_term/d_term, m_alpha), 1.0);
        assert(w_term < w_upper_bound);
        assert(d_term >= dist_lower_bound);

        if(dist(gen) < connection_prob/max_connection_prob) {
            m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
        }
    }
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::weightLayerTargetLevel(int layer) const {
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


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::partitioningBaseLevel(int layer1, int layer2) const {

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct comparison level
        assert(0 <= layer1 && layer1 < m_layers);
        assert(0 <= layer2 && layer2 < m_layers);
        auto volume_requested  = m_w0*std::pow(2,layer1+1) * m_w0*std::pow(2,layer2+1) / m_W; // v(i,j) = wi*wj/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


} // namespace girgs
