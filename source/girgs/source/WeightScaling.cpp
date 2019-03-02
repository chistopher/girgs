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
    // Conceptionally, prefix_sum contains a prefix sum of the weights vector sorted decreasingly.
    // To speed up sorting and prefix computation, we only sort as many elements as we need, hence
    // the vector contains three segments:
    //  - front element is always 0 and functions as a sentinel value to avoid boundary checks in
    //    function original_value
    //  - the prefix sum section [ps_begin, ps_end)
    //  - unordered section [ps_end, end) storing remaining elements
    std::vector<double> prefix_sum(weights.size() + 1);
    prefix_sum[0] = 0;
    const auto ps_begin = prefix_sum.begin() + 1;
          auto ps_end   = ps_begin;
    const auto end      = prefix_sum.end();
    using iter_t = decltype(ps_begin);
    auto original_value = [] (iter_t x) {return *x - *(x - 1);};

    // compute some constant stuff
    const auto n = weights.size();
    auto max_weight = 0.0;
    auto W = 0.0, sq_W = 0.0;
    {
        #pragma omp parallel for reduction(+:W, sq_W), reduction(max: max_weight)
        for (int i = 0; i < n; ++i) {
            const auto each = weights[i];
            prefix_sum[i+1] = each; // copy to prefix_sum

            W += each;
            sq_W += each * each;
            max_weight = std::max(max_weight, each);
        }
    }

    auto sorted_upper = std::numeric_limits<double>::min();
    auto prefixsum_upto = [&](double c) {
        const auto thresh = W / std::pow(2.0 * c, dimension) / max_weight;

        if (c > sorted_upper) {
            // we need more points, make sure we have more!
            assert(ps_end != prefix_sum.end());

            // move all not yet sorted elements larger than the threshold directly next to the sorted ones
            auto new_ps_end = std::partition(ps_end, end, [=](double x) { return x >= thresh; });
            if (new_ps_end == ps_end) return ps_end;

            std::sort(ps_end, new_ps_end, std::greater<double>());

            // compute the inclusive prefix sum starting one early to "connect" it to the already sorted
            // sequence. Initially we rewrite the sentinel value, so no additional checks are necessary
            std::partial_sum(ps_end - 1, new_ps_end, ps_end -1);

#ifndef NDEBUG
            for(auto it = ps_begin + 1; it != new_ps_end; ++it)
                assert(original_value(it) <= original_value(it - 1));
#endif

            ps_end = new_ps_end;
            sorted_upper = (ps_end == prefix_sum.end()) ? std::numeric_limits<double>::max() : c;
            return ps_end;

        } else {
            // the splitter is within our already sorted segment
            return std::lower_bound(ps_begin, ps_end, thresh,
                [&] (const double& x, const double thresh) {return x - *(&x - 1) > thresh;});

        }
    };

    // my function to do the exponential search on
    auto f = [ps_begin, dimension, W, sq_W, n, &prefixsum_upto, &original_value](double c) {
        // compute rich club
        const auto richclub_end = prefixsum_upto(c);

        // compute overestimation
        const auto pow2c = pow(2 * c, dimension);

        // subtract error
        auto error = 0.0;
        auto y = richclub_end - 1;
        for(auto x = ps_begin; x != richclub_end; ++x) {
            const auto fac = original_value(x) / W;
            const auto my_thres = 1.0 / fac / pow2c;

            // search smallest element larger than my_thresh; since `my_thresh` is non-decreasing,
            // y is non-decreasing either and we can use the old values as a lower-bound for the
            // next one
            for(; y >= ps_begin && original_value(y) < my_thres; y--);
            if (y < ps_begin) break;

            /**
              * sum_{k < j, k != i}{ std::pow(2*c,dimension)*(w1*w_k/W)-1.0 }
              * = sum_{k < j, k != i}{ std::pow(2*c,dimension)*(w1*w_k/W) } - j
              * = sum_{k < j}{ std::pow(2*c,dimension)*(w1*w_k/W) } - j - (0 if j < i else std::pow(2*c,dimension)*(w1*w_i/W) - 1)
              * = pow2c * w1/W * (sum_j{w_j}) - j - (0 if j < i else pow2c*(w1/W)*w_i - 1)
              */

            error += *y * pow2c * fac - (std::distance(ps_begin, y) + 1);

            if (y >= x) {
                // we have to subtract the self-contribution of x == y
                error -= pow2c * fac * original_value(x) - 1.0;
            }
        }

        const auto overestimation = pow2c * (W - sq_W / W);
        const auto res = (overestimation - error) / n;

        return res;
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