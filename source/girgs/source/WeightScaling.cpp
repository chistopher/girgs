#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <cmath>
#include <vector>

#include <girgs/WeightScaling.h>

namespace girgs {

// helper for scale weights
static double exponentialSearch(const std::function<double(double)> &f, double desiredValue, double accuracy = 0.02, double lower = 1.0,
                                double upper = 2.0) {

    // scale interval up if necessary
    while (f(upper) < desiredValue) {
        lower = upper;
        upper *= 2;
    }

    // scale interval down if necessary
    while (f(lower) > desiredValue) {
        upper = lower;
        lower /= 2;
    }

    // do binary search
    auto mid = f((upper + lower) / 2);
    while (std::abs(mid - desiredValue) > accuracy) {
        if (mid < desiredValue)
            lower = (upper + lower) / 2;
        else
            upper = (upper + lower) / 2;
        mid = f((upper + lower) / 2);
    }

    return (upper + lower) / 2;
}

// helper for scale weights
double estimateWeightScalingThreshold(const std::vector<double> &weights, double desiredAvgDegree, int dimension) {
    // compute some constant stuff
    auto max_weight = 0.0;
    auto W = 0.0, sq_W = 0.0;
    {
        const auto n = weights.size();

        #pragma omp parallel for reduction(+:W, sq_W), reduction(max: max_weight)
        for (int i = 0; i < n; ++i) {
            const auto each = weights[i];
            W += each;
            sq_W += each * each;
            max_weight = std::max(max_weight, each);
        }
    }

    // rather than sorting all weights (which can dominate the runtime for large n and small avgDeg)
    // we use stable_partition to find the "rich_club" and then only sort it. If the rich_club grows
    // over the course of the algorithm (i.e., c is growing), we apply this idea recursively on the
    // elements not sorted yet.
    std::vector<double> sorted_weights = weights;
    std::vector<double> prefix_sum(weights.size());
    auto sorted_end = sorted_weights.begin(); // points to the first element not yet sorted
    auto sorted_upper = std::numeric_limits<double>::min();
    auto get_sorted_end = [&](double c) {
        const auto thresh = W / std::pow(2.0 * c, dimension) / max_weight;

        if (c > sorted_upper) {
            assert(sorted_end != sorted_weights.end());

            // we need more points
            auto it = std::stable_partition(sorted_end, sorted_weights.end(), [=](double x) { return x >= thresh; });
            std::sort(sorted_end, it, std::greater<double>());

            // update prefixsum
            {
                // we are adding prefix (sum up to the element summed in the last iteration, if any)
                // to the first element of our new range to continue the prefix sum
                const auto prefix = (sorted_end == sorted_weights.begin()) ? 0 : *(sorted_end - 1);

                double tmp = *sorted_end + prefix;
                std::swap(*sorted_end, tmp);
                std::partial_sum(sorted_end, it, prefix_sum.begin() + std::distance(sorted_weights.begin(), sorted_end));
                std::swap(*sorted_end, tmp);
            }

            // if we reached the end of the input, we set the sorted_upper to max() in order to every
            // going into the partition branch again.
            sorted_upper = (sorted_end == sorted_weights.end()) ? std::numeric_limits<double>::max() : c;
            sorted_end = it;
            return sorted_end;

        } else {
            // the splitter is within our already sorted segment
            return std::lower_bound(sorted_weights.begin(), sorted_end, thresh, std::greater<double>());
        }
    };

    // my function to do the exponential search on
    auto f = [&](double c) {
        // compute rich club
        const auto it = get_sorted_end(c);
        const auto num_richclub = static_cast<int>(std::distance(sorted_weights.begin(), it));

        // compute overestimation
        const auto pow2c = pow(2 * c, dimension);
        const auto thresh = 1.0 / pow2c;

        // subtract error
        auto error = 0.0;
        int j = num_richclub;
        for (int i = 0; i < num_richclub; ++i) {
            const auto fac = sorted_weights[i] / W;
            const auto my_thres = thresh / fac;

            while (--j >= 0 && sorted_weights[j] <= my_thres);

            if (j < 0) break;

            const auto self_contribution = (i > j) ? 0.0 : (prefix_sum[i] * fac - 1.0);

            error += pow2c * fac * prefix_sum[j] - j - self_contribution;
        }

        const auto overestimation = pow2c * (W - sq_W / W);
        return (overestimation - error);
    };

    // do exponential search on expected average degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree * weights.size(), 0.02 * weights.size());

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
    assert(alpha != 1.0); // somehow breaks for alpha 1.0

    // compute some constant stuff
    auto W = std::accumulate(weights.begin(), weights.end(), 0.0);
    auto sum_sq_w = 0.0; // sum_{v\in V} (w_v^2/W)
    auto sum_w_a = 0.0; // sum_{v\in V} (w_v  /W)^\alpha
    auto sum_sq_w_a = 0.0; // sum_{v\in V} (w_v^2/W)^\alpha

    //   sum_{u\in V} sum_{v\in V} (wu*wv/W)^\alpha
    // = sum_{u\in V} sum_{v\in V} wu^\alpha * (wv/W)^\alpha
    // = sum_{u\in V} wu^\alpha sum_{v\in V} (wv/W)^\alpha
    auto sum_wwW_a = 0.0;
    auto max_w = 0.;

    // this loop causes >= 70% of runtime
    {
        const auto n = static_cast<int>(weights.size());
        #pragma omp parallel for reduction(+:sum_sq_w, sum_w_a, sum_sq_w_a, sum_wwW_a), reduction(max:max_w)
        for (int i = 0; i < n; ++i) {
            const auto each = weights[i];

            const auto each_W = each / W;
            const auto pow_each = pow(each, alpha);
            const auto pow_each_W = pow(each / W, alpha);

            sum_sq_w += each * each_W;
            sum_wwW_a += pow_each;
            sum_w_a += pow_each_W;
            sum_sq_w_a += pow_each * pow_each_W;
            max_w = std::max(each, max_w);
        }
    }
    sum_wwW_a *= sum_w_a;
    const auto max_w_W = max_w / W;

    const auto factor1 = (W - sum_sq_w) * (1 + 1 / (alpha - 1)) * (1 << dimension);
    const auto factor2 = pow(2, alpha * dimension) / (alpha - 1) * (sum_wwW_a - sum_sq_w_a);

    // my function to do the exponential search on
    std::vector<double> sorted_weights;
    const auto upper = 2.0;
    {
        // derivation of thresh, see below in exp-search callback
        const auto thresh = exp(dimension * log(0.5 / pow(upper, 1.0 / alpha / dimension)) - log(max_w_W));
        std::copy_if(weights.cbegin(), weights.cend(), std::back_inserter(sorted_weights), [thresh](double x) { return x > thresh; });
        sorted_weights.push_back(std::numeric_limits<double>::min()); // sentinel
        std::sort(sorted_weights.begin(), sorted_weights.end(), std::greater<double>());
    }

    std::vector<double> rich_club;
    // it's not really necessary to optimize this callback; less than 10% of computation time are spent here
    auto f = [alpha, dimension, W, factor1, factor2, max_w_W, &sorted_weights, &rich_club, upper](double c) {
        assert(c <= upper);

        auto d = dimension;
        auto a = alpha;

        // as originally in Marianne's thesis
        auto long_and_short_with_error = pow(c, 1 / a) * factor1 - c * factor2;

        // get error for long and short edges
        auto n = sorted_weights.size();
        auto short_error = 0.0;
        auto long_error = 0.0;

        // get rich club
        std::vector<double> rich_club;
        auto w_n = sorted_weights.front();

        /* We re-write the "rich-condition", s.t. no pows have to be
         * evalutated. Let w := sorted_weights[i], then:
         *
         *     pow(c, 1/a/d) * pow(w * max_w_W, 1.0/d) > 0.5
         * <=> pow(w * max_w_W, 1.0/d) > 0.5 / pow(c, 1/a/d)
         * <=> log(w * max_w_W) > d * log(0.5 / pow(c, 1/a/d))
         * <=> log(w) > d * log(0.5 / pow(c, 1/a/d)) - log(max_w_W)
         * <=> w > exp(d * log(0.5 / pow(c, 1/a/d)) - log(max_w_W))
         */
        const auto thresh = exp(d * log(0.5 / pow(c, 1 / a / d)) - log(max_w_W));
        for (int i = 0; /* break using sentinel */; ++i) {
            if (sorted_weights[i] < thresh)
                break;

            rich_club.push_back(sorted_weights[i]);
        }

        // compute errors
        const auto fac1 = (1 << d) * pow(c, 1 / a);
        const auto base = pow(0.5, d - a * d);
        for (int i = 0; i < rich_club.size(); ++i) {
            for (int j = 0; j < rich_club.size(); ++j) {
                if (i == j) continue;

                auto w_term = rich_club[i] * rich_club[j] / W;
                auto crazy_w = pow(c, 1 / a / d) * pow(w_term, 1.0 / d);
                if (crazy_w <= 0.5)
                    break;

                short_error += fac1 * w_term - 1.0;
                long_error += pow(w_term, a) * (base - pow(crazy_w, d - a * d));
            }
        }

        long_error *= c * d * (1 << d) / (d - a * d);

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
    return pow(estimated_c, 1 / alpha); // return scaling
}

} // namespace girgs