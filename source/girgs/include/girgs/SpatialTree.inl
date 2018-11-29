
#include <chrono>
#include <omp.h>

namespace girgs {


template<unsigned int D>
void SpatialTree<D>::generateEdges(std::vector<Node>& graph, double alpha, int seed) {

    // init member
    m_alpha = alpha;
    m_1 = 0;
    m_2 = 0;
    m_gen = std::mt19937(seed >= 0 ?  seed : std::random_device()());
    m_dist = std::uniform_real_distribution<>(0.0, 1.0);

    // determine min max and sum of weights
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

    // sample all edges
    // visitRoot_parallel();

    const auto num_threads = omp_get_max_threads();
    std::cout << "threads " << num_threads << '\n';
	if (num_threads == 1) {
		auto start1 = std::chrono::high_resolution_clock::now();
		visitCellPair(0, 0, 0);
		auto start2 = std::chrono::high_resolution_clock::now();
		std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start2 - start1).count() << '\n';
		return;
	}
    const unsigned int first_parallel_level = std::ceil(std::log2(4.0*num_threads) / D);
    const auto parallel_cells = SpatialTreeCoordinateHelper<D>::numCellsInLevel(first_parallel_level);
    const auto first_parallel_cell = SpatialTreeCoordinateHelper<D>::firstCellOfLevel(first_parallel_level);
    assert(first_parallel_level < m_levels);

	auto start1 = std::chrono::high_resolution_clock::now();
    auto parallel_calls = std::vector<std::vector<unsigned int>>(parallel_cells);
    visitCellPair_sequentialStart(0,0,0, first_parallel_level, parallel_calls);
	auto start2 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic), num_threads(num_threads)
    for(int i=0; i< parallel_cells; ++i) {
        auto current_cell = first_parallel_cell + i;
        for(auto each : parallel_calls[i])
            visitCellPair(current_cell, each, first_parallel_level);
    }
	auto start3 = std::chrono::high_resolution_clock::now();

	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start2 - start1).count() << '\n';
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start3 - start2).count() << '\n';

    // after sampling the graph m_1+m_2 should always be n(n-1) to ensure that all edges were considered
    assert(m_1 + m_2 == graph.size()*(graph.size()-1));
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
void SpatialTree<D>::visitRoot_parallel() {
    using Helper = SpatialTreeCoordinateHelper<D>;

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[0]){
        assert(partitioningBaseLevel(layer_pair.first, layer_pair.second) == 0);
        if(layer_pair.first <= layer_pair.second)
            sampleTypeI(0, 0, 0, layer_pair.first, layer_pair.second);
    }

    // recursive call for all children pairs (a,b) of root node 0
    #pragma omp parallel for schedule(dynamic), num_threads(3), if(D>1)
    for(int a = Helper::firstChild(0); a<=Helper::lastChild(0); ++a)
        for (auto b = a; b <= Helper::lastChild(0); ++b)
            visitCellPair(a, b, 1);
}


template<unsigned int D>
void SpatialTree<D>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {
    using Helper = SpatialTreeCoordinateHelper<D>;

    auto touching = m_helper.touching(cellA, cellB, level);

    // TODO consider early return if A or B empty

    if(cellA == cellB || touching) {

        // sample all type 1 occurrences with this cell pair
        for(auto& layer_pair : m_layer_pairs[level]){
            assert(partitioningBaseLevel(layer_pair.first, layer_pair.second) == level);
            if(cellA != cellB || layer_pair.first <= layer_pair.second)
                sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
        }

    } else { // not touching

        // sample all type 2 occurrences with this cell pair
        // TODO consider early break 1st for loop
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
void SpatialTree<D>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{

    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);

    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =(cellA == cellB && i==j ? kA+1 : 0); kB<sizeV_j_B; ++kB) {
            Node* nodeInA = m_weight_layers[i].kthPoint(cellA, level, kA);
            Node* nodeInB = m_weight_layers[j].kthPoint(cellB, level, kB);

            // points are in correct cells
            assert(cellA == m_helper.cellForPoint(nodeInA->coord, level));
            assert(cellB == m_helper.cellForPoint(nodeInB->coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA->weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB->weight/m_w0)));

            assert(nodeInA->index != nodeInB->index);
            //#pragma omp atomic update
            //m_1 += 1 + (nodeInA->index != nodeInB->index);
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
    auto sizeV_i_A = m_weight_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_weight_layers[j].pointsInCell(cellB, level);
    if(sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

    //#pragma omp atomic update
    //m_2 += 2*sizeV_i_A*sizeV_j_B;

    // implicit sampling
    auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    auto dist_lower_bound = std::pow(m_helper.dist(cellA, cellB, level), dimension);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    if(max_connection_prob <= 1e-10)
        return;
    auto geo = [this](double p) -> long long {
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
        auto connection_prob = std::min(std::pow(w/d, m_alpha), 1.0);
        if(m_dist(m_gen) < connection_prob/max_connection_prob) {
            nodeInA->edges.push_back(nodeInB);
            //nodeInB->edges.push_back(nodeInA);
        }
        r += geo(max_connection_prob);
    }
}


template<unsigned int D>
unsigned int SpatialTree<D>::weightLayerTargetLevel(int layer) const {
    // -1 coz w0 is the upper bound for layer 0 in paper and our layers are shifted by -1
    auto result = std::max((m_baseLevelConstant - layer - 1) / (int)D, 0);
    {   // a lot of assertions that we have the correct insertion level
        assert(0 <= layer && layer < m_layers);
        assert(0 <= result && result <= m_levels); // note the result may be one larger than the deepest level (hence the <= m_levels)
        auto volume_requested  = m_w0*m_w0*std::pow(2,layer+1)/m_W; // v(i) = w0*wi/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
    return static_cast<unsigned int>(result);
}


template<unsigned int D>
unsigned int SpatialTree<D>::partitioningBaseLevel(int layer1, int layer2) const {
    /*
     * v(i,j) -- partitioning base level
     *
     * w0 is min weight
     * wi = 2*w(i-1) = 2^i * w0 // split wi into power of two and w0
     *
     * 2^(-ld) = wi*wj/W
     * = 2^i*w0 * 2^j*w0 / W
     * = 2^i * 2^j / (W/w0^2)
     * -> -ld = i+j - log2(W/w0^2)
     * -> l  = (log2(W/w0^2) - (i+j)) / d
     *
     * a point with weight w is inserted into layer floor(log2(w/w0))
     *  -> w0 is inserted in layer 0 (instead of layer 1 like in paper)
     *  -> so in fact wi in our implementation equals w(i+1) in paper
     * a pair of layers is compared in level (log2(W/w0^2) - (i+j+2)) / d rounded down
     * - +2 to shift from our wi back to paper wi
     * - rounding down means a level with less depth like requested in paper
    */

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
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
    return static_cast<unsigned int>(result);
}


template<unsigned int D>
bool SpatialTree<D>::checkEdgeExplicit(double dist, double w1, double w2) {
    auto w_term = w1*w2/m_W;
    auto d_term = std::pow(dist, D);
    if(m_alpha == std::numeric_limits<double>::infinity())
        return d_term < w_term;

    auto edge_prob = std::min(std::pow(w_term/d_term, m_alpha), 1.0);
    return m_dist(m_gen) < edge_prob;
}



} // namespace girgs
