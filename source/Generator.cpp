
#include <Generator.h>

#include <fstream>
#include <iomanip>
#include <random>

#include <SpatialTree.h>



void Generator::setWeights(const std::vector<double>& weights) {
    auto n = weights.size();
    assert(m_graph.empty() || m_graph.size() == n);
    if(m_graph.empty()) m_graph.resize(n);
    for(int i=0; i<n; ++i)
        m_graph[i].weight = weights[i];
}


void Generator::setWeights(int n, double ple, int weightSeed) {
    assert(m_graph.empty() || m_graph.size() == n);
    assert(ple <= -2);
    if(m_graph.empty()) m_graph.resize(n);

    auto gen = std::mt19937(weightSeed >= 0 ?  weightSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        m_graph[i].weight = std::pow(std::pow(n,ple+1)*dist(gen) + 1, 1.0/(ple+1));
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
    alpha = std::numeric_limits<double>::infinity(); // TODO find a way to actually consider alpha

    { // TODO change the whole block
        // shamefully copy all weights
        auto currentWeights = std::vector<double>(n);
        for(int i=0; i<n; ++i)
            currentWeights[i] = m_graph[i].weight;
        // scale weights
        auto scaling = estimateWeightScaling(currentWeights, desiredAvgDegree, dimension, alpha);
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


void Generator::generateTreshold() {
    generate(std::numeric_limits<double>::infinity(), 0);
}


const std::vector<Node>& Generator::generate(
        int n, int dimension, double ple, double alpha, int desiredAvgDegree, int weightSeed, int positionSeed, int samplingSeed) {

    setWeights(n, ple, weightSeed);
    setPositions(n,dimension, positionSeed);

    scaleWeights(desiredAvgDegree, dimension, alpha);

    generate(alpha, samplingSeed);

    return m_graph;
}



double Generator::avg_degree() const {
    auto edges = 0.0;
    for(auto& each : graph())
        edges += each.edges.size();
    return edges / m_graph.size();
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
            if(each.index < neighbor->index)
                f << '\t' << each.index << "\t-- " << neighbor->index << ";\n";
    }
    f << "}\n";
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


double Generator::estimateWeightScaling(const std::vector<double>& weights, int desiredAvgDegree, int dimension, double alpha) const {

    using namespace std;

    // currently only for threshold model!!!!
    assert(alpha == std::numeric_limits<double>::infinity());

    // compute some constant stuff
    auto max_weight = *max_element(weights.begin(), weights.end());
    auto W = 0.0, sq_W = 0.0;
    for(auto each : weights){
        W += each;
        sq_W += each*each;
    }

    // my function to do the exponential search on
    auto f = [W, sq_W, &weights, dimension, n(weights.size()), max_weight](double c) {
        // compute rich club
        vector<double> rich_club;
        for(auto weight : weights)
            if(std::pow(2*c,dimension) * (weight*max_weight/W) > 1.0)
                rich_club.push_back(weight);
        sort(rich_club.begin(), rich_club.end(), greater<>());
        // compute overestimation
        auto overestimation = pow(2, dimension) * pow(c, dimension) / n * (W - sq_W/W);
        // subtract error
        auto error = 0.0;
        for(int i = 0; i<rich_club.size(); ++i)
            for(int j = 0; j<rich_club.size(); ++j) {
                if(i==j) continue;
                auto w1 = rich_club[i];
                auto w2 = rich_club[j];
                auto e = max( std::pow(2*c,dimension)*(w1*w2/W)-1.0, 0.0);
                error += e;
                if(e <= 0) break;
            }
        return overestimation - error/n;
    };

    // do exponential search on avg_degree function
    auto upper = pow( weights.size()*desiredAvgDegree / pow(2.0,dimension) / (W-sq_W/W), 1.0/dimension);
    auto lower = upper / 2.0;

    // scale interval up if necessary
    while(f(upper) < desiredAvgDegree){
        lower = upper;
        upper *= 2;
    }

    // scale interval down if necessary
    while(f(lower) > desiredAvgDegree){
        upper = lower;
        lower /= 2;
    }

    // do binary search
    auto mid = f((upper+lower)/2);
    while(abs(mid - desiredAvgDegree) > 0.02) {
        if(mid < desiredAvgDegree)
            lower = (upper+lower)/2;
        else
            upper = (upper+lower)/2;
        mid = f((upper+lower)/2);
    }

    auto estimated_c = (upper+lower)/2;
    /*
     * edge iff dist < c(wi*wj/W)^(1/d)
     *
     * c(wi*wj/W)^(1/d)
     * = ( c^d * wi*wj/W )^(1/d)
     * = ( c^d*wi * c^d*wi / (c^d*W) )^(1/d)
     *
     * so we can just scale all weights by c^d
     *
     * TODO non threshold case
     */

    return pow(estimated_c, dimension); // return scaling
}
