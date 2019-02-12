#pragma once

#include <array>
#include <cmath>
#include <limits>


namespace hypergirgs {

/**
 * The filter should be used to compute bounds on the hyperbolic distance (or rather the cosh of it),
 * such that vertices with this distance are guaranteed to be connected/disconnected.
 *
 * @tparam stages
 *  The resolution of the filter (i.e. number of entries).
 */
template<size_t stages>
class DistanceFilter {
public:
    /**
     * Creates a distance filter for the hyperbolic connection probability.
     *
     * @param max_prob
     *  If max_prob < 1.0, then queries to the filter are restricted to arguments in [0..max_prob).
     *  This increases the resolution in the specified range.
     * @param R
     *  The radius of the hyperbolic disc.
     * @param T
     *  The temperature of the graph.
     */
    DistanceFilter(double max_prob, double R, double T)
    : max_connection_prob(max_prob)
    , stage_width(stages/max_prob)
    {
        for(int i=1; i < stages; i++)
            filter_stages[i] = cosh(invConnectionProb(max_connection_prob / stages * i, R, T));
        // TODO maybe we should not save these values
        filter_stages[0] = std::numeric_limits<double>::infinity(); // at distance inf the conn prob is 0; cosh(inf)=inf
        filter_stages[stages] = 1.0; // at distance 0 the conn prob is largest; cosh(0)=1
    }

    /// upper bound for cosh(dist) to have an edge probability lower than prob
    /// a vertex pair is surely disconnected if cosh(dist) is larger than the returned value
    double coshDistForProb_upperBound(double prob) const {
        assert(prob < max_connection_prob);
        const auto index = static_cast<int>(stage_width * prob);
        return filter_stages[index];
    }

    /// lower bound for cosh(dist) to have an edge probability larger than prob
    /// a vertex pair is surely connected if cosh(dist) is smaller than the returned value
    double coshDistForProb_lowerBound(double prob) const {
        assert(prob < max_connection_prob);
        const auto index = static_cast<int>(std::ceil(stage_width * prob));
        return filter_stages[index];
    }


    const double max_connection_prob;     ///< range of the filter. See c'tor
private:

    /// inverse of the edge probability function
    /// invConnectionProb(p_uv) = dist_uv
    double invConnectionProb(double p, double R, double T) const {
        return R + 2*T*std::log(1.0 / p - 1);
    }

    std::array<double, stages+1> filter_stages; ///< pos i holds the cosh(distance) at which to points would have a connection prob of (i/stages)
    const double stage_width;
};


} // namespace hypergirgs
