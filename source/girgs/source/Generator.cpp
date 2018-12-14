
#include <girgs/Generator.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>

#include <girgs/SpatialTree.h>


using namespace girgs;


void Generator::setWeights(const std::vector<double>& weights) {
    auto n = weights.size();
    assert(m_graph.empty() || m_graph.size() == n);
    if(m_graph.empty()) m_graph.resize(n);
    for(auto i=0u; i<n; ++i)
        m_graph[i].weight = weights[i];
}


void Generator::setWeights(int n, double ple, int weightSeed) {
    assert(m_graph.empty() || m_graph.size() == n);
    assert(ple <= -2);
    if(m_graph.empty()) m_graph.resize(n);

    auto gen = std::mt19937(weightSeed >= 0 ?  weightSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        m_graph[i].weight = std::pow((std::pow(n,ple+1)-1)*dist(gen) + 1, 1/(ple+1));
}


void Generator::setPositions(const std::vector<std::vector<double>>& positions) {
    auto n = positions.size();
    assert(m_graph.empty() || m_graph.size() == n);
    if(m_graph.empty()) m_graph.resize(n);
    for(int i=0; i<n; ++i) {
        assert(positions[i].size() == positions.front().size()); // all same dimension
        m_graph[i].coord = positions[i];
    }
}


void Generator::setPositions(int n, int dimension, int positionSeed) {
    assert(m_graph.empty() || m_graph.size() == n);
    if(m_graph.empty()) m_graph.resize(n);

    auto gen = std::mt19937(positionSeed >= 0 ?  positionSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)

    for(int i=0; i<n; ++i) {
        m_graph[i].coord.resize(dimension);
        for (int d=0; d<dimension; ++d)
            m_graph[i].coord[d] = dist(gen);
    }
}


double Generator::scaleWeights(int desiredAvgDegree, int dimension, double alpha) {
    assert(!m_graph.empty());
    auto n = m_graph.size();

    { // TODO change the whole block
        // shamefully copy all weights
        auto currentWeights = std::vector<double>(n);
        for(int i=0; i<n; ++i)
            currentWeights[i] = m_graph[i].weight;

        // estimate scaling with binary search
        double scaling;
        if(alpha > 10.0)
            scaling = estimateWeightScalingThreshold(currentWeights, desiredAvgDegree, dimension);
        else if(alpha > 0.0 && alpha != 1.0)
            scaling = estimateWeightScaling(currentWeights, desiredAvgDegree, dimension, alpha);
        else
            throw("I do not know how to scale weights for desired alpha :(");

        // scale weights
        for(auto& each : m_graph)
            each.weight *= scaling;
        return scaling;
    }
}




void Generator::generate(double alpha, int samplingSeed) {
    assert(!m_graph.empty());
    for(auto& each : m_graph) each.edges.clear();
    auto dimension = m_graph.front().coord.size();
    switch(dimension) {
        case 1: SpatialTree<1>().generateEdges(m_graph, alpha, samplingSeed); break;
        case 2: SpatialTree<2>().generateEdges(m_graph, alpha, samplingSeed); break;
        case 3: SpatialTree<3>().generateEdges(m_graph, alpha, samplingSeed); break;
        case 4: SpatialTree<4>().generateEdges(m_graph, alpha, samplingSeed); break;
        case 5: SpatialTree<5>().generateEdges(m_graph, alpha, samplingSeed); break;
        default:
            std::cout << "Dimension " << dimension << " not supported." << std::endl;
            std::cout << "No edges generated." << std::endl;
            break;
    }
}


std::vector<Node> Generator::generate(
        int n, int dimension, double ple, double alpha, int desiredAvgDegree, int weightSeed, int positionSeed, int samplingSeed) {

    setWeights(n, ple, weightSeed);
    setPositions(n,dimension, positionSeed);

    scaleWeights(desiredAvgDegree, dimension, alpha);

    generate(alpha, samplingSeed);

    return std::move(m_graph);
}


void Generator::generateThreshold() {
    generate(std::numeric_limits<double>::infinity(), 0);
}


double Generator::avg_degree() const {
    return 2.0 * edges() / m_graph.size();
}

unsigned int girgs::Generator::edges() const {
    auto edges = 0u;
    for (auto& each : graph())
        edges += each.edges.size();
    return edges;
}


void Generator::saveDot(std::string file) const {
    if(m_graph.empty()){
        std::cout << "no graph generated" << std::endl;
        return;
    }

    auto f = std::ofstream(file);
    f << "graph girg {\n\toverlap=scale;\n\n";
    for(auto& each : m_graph){
        f   << '\t' << each.index << " [label=\""
            << std::setprecision(2) << std::fixed << each.weight << std::defaultfloat << std::setprecision(6)
            << "\", pos=\"";
        for(auto d=0u; d<each.coord.size(); ++d)
            f << (d==0 ? "" : ",") << each.coord[d];
        f << "\"];\n";
    }
    f << '\n';
    for(auto& each : m_graph){
        for(auto neighbor : each.edges)
            f << '\t' << each.index << "\t-- " << neighbor->index << ";\n";
    }
    f << "}\n";
}


void girgs::Generator::saveEdgeList(std::string file) const {
    auto f = std::ofstream(file);
    f << m_graph.size() << ' ' << edges() << '\n';
    for (auto& from : m_graph)
        for (auto to : from.edges)
            f << from.index << ' ' << to->index << '\n';
}



std::vector<double> Generator::weights() const {
    auto result = std::vector<double>(m_graph.size());
    for(int i=0; i<m_graph.size(); ++i)
        result[i] = m_graph[i].weight;
    return result;
}

std::vector<std::vector<double>> Generator::positions() const {
    auto result = std::vector<std::vector<double>>(m_graph.size());
    for(int i=0; i<m_graph.size(); ++i)
        result[i] = m_graph[i].coord;
    return result;
}


double Generator::estimateWeightScalingThreshold(const std::vector<double>& weights, int desiredAvgDegree, int dimension) const {

    // compute some constant stuff
    auto max_weight = *std::max_element(weights.begin(), weights.end());
    auto W = 0.0, sq_W = 0.0;
    for(auto each : weights){
        W += each;
        sq_W += each*each;
    }

    // my function to do the exponential search on
    auto f = [W, sq_W, &weights, dimension, max_weight](double c) {
        // compute rich club
        std::vector<double> rich_club;
        for(auto weight : weights)
            if(std::pow(2*c,dimension) * (weight*max_weight/W) > 1.0)
                rich_club.push_back(weight);
        sort(rich_club.begin(), rich_club.end(), std::greater<double>());
        // compute overestimation
        auto overestimation = pow(2, dimension) * pow(c, dimension) * (W - sq_W/W);
        // subtract error
        auto error = 0.0;
        for(int i = 0; i<rich_club.size(); ++i)
            for(int j = 0; j<rich_club.size(); ++j) {
                if(i==j) continue;
                auto w1 = rich_club[i];
                auto w2 = rich_club[j];
                auto e = std::max( std::pow(2*c,dimension)*(w1*w2/W)-1.0, 0.0);
                error += e;
                if(e <= 0) break;
            }
        return (overestimation - error) / weights.size();
    };

    // do exponential search on expected average degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree);

    /*
     * edge iff dist < c(wi*wj/W)^(1/d)
     *
     * c(wi*wj/W)^(1/d)
     * = ( c^d * wi*wj/W )^(1/d)
     * = ( c^d*wi * c^d*wi / (c^d*W) )^(1/d)
     *
     * so we can just scale all weights by c^d
     */
    return pow(estimated_c, dimension); // return scaling
}

double Generator::estimateWeightScaling(const std::vector<double> &weights, int desiredAvgDegree, int dimension, double alpha) const {

    using namespace std;

    assert(alpha != 1.0); // somehow breaks for alpha 1.0

    // compute some constant stuff
    auto W = std::accumulate(weights.begin(), weights.end(), 0.0);
    auto sum_sq_w   = 0.0; // sum_{v\in V} (w_v^2/W)
    auto sum_w_a    = 0.0; // sum_{v\in V} (w_v  /W)^\alpha
    auto sum_sq_w_a = 0.0; // sum_{v\in V} (w_v^2/W)^\alpha
    for(auto each : weights){
        sum_sq_w   += each*each/W;
        sum_w_a    += pow(each/W, alpha);
        sum_sq_w_a += pow(each*each/W, alpha);
    }

    //   sum_{u\in V} sum_{v\in V} (wu*wv/W)^\alpha
    // = sum_{u\in V} sum_{v\in V} wu^\alpha * (wv/W)^\alpha
    // = sum_{u\in V} wu^\alpha sum_{v\in V} (wv/W)^\alpha
    auto sum_wwW_a = 0.0;
    for(auto each : weights)
        sum_wwW_a += pow(each, alpha)*sum_w_a;

    auto factor1 = (W-sum_sq_w) * (1+1/(alpha-1)) * (1<<dimension);
    auto factor2 = pow(2, alpha*dimension) / (alpha-1) * (sum_wwW_a - sum_sq_w_a);

    // my function to do the exponential search on
    auto sorted_weights = weights;
    sort(sorted_weights.begin(), sorted_weights.end(), greater<double>());

    auto f = [alpha, dimension, W, factor1, factor2, &sorted_weights](double c) {
        auto d = dimension;
        auto a = alpha;

        // as originally in Marianne's thesis
        auto long_and_short_with_error = pow(c, 1/a) * factor1 - c * factor2;

        // get error for long and short edges
        auto n = sorted_weights.size();
        auto short_error = 0.0;
        auto long_error = 0.0;

        // get rich club
        vector<double> rich_club;
        auto w_n = sorted_weights.front();
        for(int i=0; i<n; ++i){
            auto crazy_w = pow(c, 1/a/d) * pow(sorted_weights[i] * w_n / W, 1.0/d);
            if(crazy_w > 0.5)
                rich_club.push_back(sorted_weights[i]);
            else
                break;
        }

        // compute errors
        for(int i = 0; i<rich_club.size(); ++i) {
            for (int j = 0; j < rich_club.size(); ++j) {
                if (i == j) continue;

                auto w_term = rich_club[i] * rich_club[j] / W;
                auto crazy_w = pow(c, 1 / a / d) * pow(w_term, 1.0 / d);
                if (crazy_w <= 0.5)
                    break;

                short_error += (1<<d)*pow(c, 1/a)*w_term -1.0;
                long_error += c * pow(w_term, a) * d * (1 << d) / (d - a * d) * (pow(0.5, d - a * d) - pow(crazy_w, d - a * d));
            }
        }
        return (long_and_short_with_error - short_error - long_error)/n;
    };

    // do exponential search on avg_degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree);

    /*
     * Pr(edge) = Pr(c * 1/dist^ad * (wi*wj/W)^a )
     *
     * c * (wi*wj/W)^a
     * = (c^{1/a} wi*wj/W)^a
     * = (c^{1/a}wi* (c^{1/a}wj / (c^{1/a}W)^a
     *
     * so we can just scale all weights by (c^{1/a}
     */

    return pow(estimated_c, 1/alpha); // return scaling
}

double Generator::exponentialSearch(std::function<double(double)> f, double desiredValue, double accuracy, double lower, double upper) const {

    // scale interval up if necessary
    while(f(upper) < desiredValue){
        lower = upper;
        upper *= 2;
    }

    // scale interval down if necessary
    while(f(lower) > desiredValue){
        upper = lower;
        lower /= 2;
    }

    // do binary search
    auto mid = f((upper+lower)/2);
    while(std::abs(mid - desiredValue) > accuracy) {
        if(mid < desiredValue)
            lower = (upper+lower)/2;
        else
            upper = (upper+lower)/2;
        mid = f((upper+lower)/2);
    }

    return (upper+lower)/2;
}
