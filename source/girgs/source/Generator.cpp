
#include <girgs/Generator.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <mutex>

#include <omp.h>

#include <girgs/SpatialTree.h>


namespace girgs {

// helper for scale weights
static double exponentialSearch(const std::function<double(double)>& f, double desiredValue, double accuracy = 0.02, double lower = 1.0, double upper = 2.0) {

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

// helper for scale weights
double estimateWeightScalingThreshold(const std::vector<double>& weights, double desiredAvgDegree, int dimension) {
    // compute some constant stuff
    auto max_weight = 0.0;
    auto W = 0.0, sq_W = 0.0;
    {
        const auto n = weights.size();

        #pragma omp parallel for reduction(+:W, sq_W), reduction(max: max_weight)
        for(int i = 0; i < n; ++i) {
            const auto each = weights[i];
            W += each;
            sq_W += each * each;
            max_weight = std::max(max_weight, each);
        }
    }

    // filter and sort candidates for rich_club
    const double upper = 2.0;
    std::vector<double> sorted_weights;
    {
        const auto thresh = W / std::pow(2*upper,dimension) / max_weight;
        std::copy_if(weights.cbegin(), weights.cend(), std::back_inserter(sorted_weights),
                     [thresh] (double x) {return x > thresh;});
        sort(sorted_weights.begin(), sorted_weights.end(), std::greater<double>());
        sorted_weights.push_back(std::numeric_limits<double>::min());
    }

    // my function to do the exponential search on
    std::vector<double> rich_club;
    auto f = [&] (double c) {
        // compute rich club
        rich_club.clear();
        const auto thresh = W / std::pow(2*c,dimension) / max_weight;

        for(int i=0; /* sentinel based */; ++i) {
            if (sorted_weights[i] < thresh)
                break;

            rich_club.push_back(sorted_weights[i]);
        }

        // compute overestimation
        auto overestimation = pow(2, dimension) * pow(c, dimension) * (W - sq_W/W);

        // subtract error
        auto error = 0.0;
        for(int i = 0; i<rich_club.size(); ++i) {
            for (int j = 0; j < rich_club.size(); ++j) {
                if (i == j)
                    continue;

                const auto w1 = rich_club[i];
                const auto w2 = rich_club[j];
                const auto e = std::max(std::pow(2 * c, dimension) * (w1 * w2 / W) - 1.0, 0.0);
                error += e;

                if (e <= 0)
                    break;
            }
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


// helper for scale weights
double estimateWeightScaling(const std::vector<double> &weights, double desiredAvgDegree, int dimension, double alpha) {

    using namespace std;

    assert(alpha != 1.0); // somehow breaks for alpha 1.0

    // compute some constant stuff
    auto W = std::accumulate(weights.begin(), weights.end(), 0.0);
    auto sum_sq_w   = 0.0; // sum_{v\in V} (w_v^2/W)
    auto sum_w_a    = 0.0; // sum_{v\in V} (w_v  /W)^\alpha
    auto sum_sq_w_a = 0.0; // sum_{v\in V} (w_v^2/W)^\alpha

    //   sum_{u\in V} sum_{v\in V} (wu*wv/W)^\alpha
    // = sum_{u\in V} sum_{v\in V} wu^\alpha * (wv/W)^\alpha
    // = sum_{u\in V} wu^\alpha sum_{v\in V} (wv/W)^\alpha
    auto sum_wwW_a = 0.0;
    auto max_w = 0.;

    {
        const auto n = static_cast<int>(weights.size());
        #pragma omp parallel for reduction(+:sum_sq_w,sum_w_a,sum_sq_w_a,sum_wwW_a), reduction(max:max_w)
        for(int i = 0; i < n; ++i) {
            const auto each = weights[i];

            const auto each_W     = each / W;
            const auto pow_each   = pow(each, alpha);
            const auto pow_each_W = pow(each / W, alpha);

            sum_sq_w   += each * each_W;
            sum_wwW_a  += pow_each;
            sum_w_a    += pow_each_W;
            sum_sq_w_a += pow_each * pow_each_W;
            max_w       = std::max(each, max_w);
        }
    }
    sum_wwW_a *= sum_w_a;
    const auto max_w_W = max_w / W;

    const auto factor1 = (W-sum_sq_w) * (1+1/(alpha-1)) * (1<<dimension);
    const auto factor2 = pow(2, alpha*dimension) / (alpha-1) * (sum_wwW_a - sum_sq_w_a);

    // my function to do the exponential search on
    std::vector<double> sorted_weights;
    const auto upper = 2.0;
    {
        // derivation of thresh, see below in exp-search callback
        const auto thresh = exp(dimension * log(0.5 / pow(upper, 1.0 / alpha / dimension)) - log(max_w_W));

        for(auto w : weights) {
            if(w > thresh)
                sorted_weights.push_back(w);
        }

        sorted_weights.push_back(std::numeric_limits<double>::min()); // sentinel
        sort(sorted_weights.begin(), sorted_weights.end(), greater<double>());
    }

    std::vector<double> rich_club;
    auto f = [alpha, dimension, W, factor1, factor2, max_w_W, &sorted_weights, &rich_club, upper](double c) {
        assert(c <= upper);

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

        /**
         * w := sorted_weights[i]
         *     pow(c, 1/a/d) * pow(w * max_w_W, 1.0/d) > 0.5
         * <=> pow(w * max_w_W, 1.0/d) > 0.5 / pow(c, 1/a/d)
         * <=> log(w * max_w_W) > d * log(0.5 / pow(c, 1/a/d))
         * <=> log(w) > d * log(0.5 / pow(c, 1/a/d)) - log(max_w_W)
         * <=> w > exp(d * log(0.5 / pow(c, 1/a/d)) - log(max_w_W))
         */
        const auto thresh = exp(d * log(0.5 / pow(c, 1/a/d)) - log(max_w_W));
        for(int i=0; /* break using sentinel */; ++i) {
            if (sorted_weights[i] < thresh)
                break;

            rich_club.push_back(sorted_weights[i]);
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
        return (long_and_short_with_error - short_error - long_error);
    };

    // do exponential search on avg_degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree * weights.size());

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


std::vector<double> generateWeights(int n, double ple, int weightSeed) {
    auto result = std::vector<double>(n);
    auto gen = std::mt19937(weightSeed >= 0 ?  weightSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        result[i] = std::pow((std::pow(n,-ple+1)-1)*dist(gen) + 1, 1/(-ple+1));
    return result;
}


std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed) {
    auto result = std::vector<std::vector<double>>(n, std::vector<double>(dimension));
    auto gen = std::mt19937(positionSeed >= 0 ?  positionSeed : std::random_device()());
    std::uniform_real_distribution<> dist; // [0..1)
    for(int i=0; i<n; ++i)
        for (int d=0; d<dimension; ++d)
            result[i][d] = dist(gen);
    return result;
}


double scaleWeights(std::vector<double>& weights, double desiredAvgDegree, int dimension, double alpha) {

    // estimate scaling with binary search
    double scaling;
    if(alpha > 10.0)
        scaling = estimateWeightScalingThreshold(weights, desiredAvgDegree, dimension);
    else if(alpha > 0.0 && alpha != 1.0)
        scaling = estimateWeightScaling(weights, desiredAvgDegree, dimension, alpha);
    else
        throw("I do not know how to scale weights for desired alpha :(");

    // scale weights
    for(auto& each : weights) each *= scaling;
    return scaling;
}

std::vector<std::pair<int, int>> generateEdges(const std::vector<double> &weights, const std::vector<std::vector<double>> &positions,
        double alpha, int samplingSeed) {

    using edge_vector = std::vector<std::pair<int, int>>;
    edge_vector result;

    std::vector<std::pair<
            edge_vector,
            uint64_t[31] /* avoid false sharing */
    > > local_edges(omp_get_max_threads());

    constexpr auto block_size = size_t{1} << 20;

    std::mutex m;
    auto flush = [&] (const edge_vector& local) {
        std::lock_guard<std::mutex> lock(m);
        result.insert(result.end(), local.cbegin(), local.cend());
    };

    auto addEdge = [&](int u, int v, int tid) {
        auto& local = local_edges[tid].first;
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
        }
    };

    auto dimension = positions.front().size();

    switch(dimension) {
        case 1: makeSpatialTree<1>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 2: makeSpatialTree<2>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 3: makeSpatialTree<3>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 4: makeSpatialTree<4>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        case 5: makeSpatialTree<5>(weights, positions, alpha, addEdge).generateEdges(samplingSeed); break;
        default:
            std::cout << "Dimension " << dimension << " not supported." << std::endl;
            std::cout << "No edges generated." << std::endl;
            break;
    }

    for(const auto& v : local_edges)
        flush(v.first);

    return result;
}


void saveDot(const std::vector<double> &weights, const std::vector<std::vector<double>> &positions,
             std::vector<std::pair<int, int>> graph, std::string file) {

    auto f = std::ofstream(file);
    f << "graph girg {\n\toverlap=scale;\n\n";
    for (int i = 0; i < weights.size(); ++i) {
        f << '\t' << i << " [label=\""
          << std::setprecision(2) << std::fixed << weights[i] << std::defaultfloat << std::setprecision(6)
          << "\", pos=\"";
        for (auto d = 0u; d < positions[i].size(); ++d)
            f << (d == 0 ? "" : ",") << positions[i][d];
        f << "\"];\n";
    }
    f << '\n';
    for (auto &edge : graph)
        f << '\t' << edge.first << "\t-- " << edge.second << ";\n";
    f << "}\n";
}

} // namespace girgs
