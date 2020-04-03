#include <algorithm>
#include <cassert>
#include <vector>

#include <gtest/gtest.h>
#include <boost/math/distributions/binomial.hpp>

#include <hypergirgs/HyperbolicTree.h>


class EdgeProbabilities_test: public ::testing::TestWithParam< std::tuple<unsigned, unsigned, double, double> > {
protected:
    const int sampleSeed = 12456;
    const int edgesSeed = 1400;

    static size_t edge2idx(unsigned u, unsigned v, size_t n) {
        const auto mm = std::minmax(u, v);
        return mm.first * n + mm.second;
    }

    static std::vector<uint64_t>
    count_edges(const std::vector<double> &radii, const std::vector<double> &angles, double T, double R, int seed,
                const unsigned num_graphs) {
        const auto n = radii.size();
        assert(radii.size() == angles.size());
        assert(n > 0);

        std::vector<uint64_t> counts(n * n, 0);
        if (!num_graphs)
            return counts;

        auto count_callback = [&counts, n](int u, int v, int /*tid*/) { counts[edge2idx(u, v, n)]++; };
        auto generator = hypergirgs::makeHyperbolicTree(radii, angles, T, R, count_callback, false);

        for (unsigned i = 0; i < num_graphs; ++i)
            generator.generate(i * 1234567 + 123);

        return counts;
    }

    static uint64_t count_violations(const std::vector<double> &radii, const std::vector<double> &angles, double T, double R,
                                     const std::vector<uint64_t> &counts, const unsigned num_graphs, const double confidence) {
        using namespace boost::math;

        const auto n = radii.size();
        assert(radii.size() == angles.size());
        assert(n > 0);

        uint64_t violations = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+:violations)
        for (int u = 0; u < n; ++u) {
            for (int v = u + 1; v < n; ++v) {
                const auto idx = edge2idx(u, v, n);
                const auto dist = hypergirgs::hyperbolicDistance(radii[u], angles[u], radii[v], angles[v]);
                const auto prob = 1.0 / (1.0 + std::exp(0.5 / T * (dist - R)));
                const auto count = counts[idx];
                const auto mean = prob * num_graphs;

                const double pvalue = 2.0 * ((count < mean) ? cdf(binomial(num_graphs, prob), count) : cdf(
                    complement(binomial(num_graphs, prob), count - 1)));

                EXPECT_GT(pvalue, confidence)
                                << "Edge (" << u << ", " << v << ") with r1=" << radii[u] << ", phi1" << angles[u] << ", r2=" << radii[v]
                                << ", phi2" << angles[v] << ", dist=" << dist << ", prob=" << prob << ", count=" << count << ", p-value="
                                << pvalue;

                violations += pvalue < confidence;
            }
        }

        return violations;
    }
};

/*
 * For positive temperatures we carry out statistical tests for each edge. Our null hypothesis is that our generator works
 * correct. Then we compute numGraphs many graph instances for the same set of points and count for each possible edge (u,v) the
 * number N(u,v) of times it is generated. Given that (u,v) has a probability of p(u,v) we expect N(u,v) to be distributed as
 * Bernoulli(numGraphs, p(u,v)). We compute its p-value and reject with a significance level of sig (corrected for the fact,
 * the we carry out (n over 2) trials -- one for each possible edge).
 */
TEST_P(EdgeProbabilities_test, StatTestCorrectProbs) {
    unsigned n, deg;
    double alpha, T;
    std::tie(n, deg, alpha, T) = GetParam();
    ASSERT_GT(T, 0.0);

    const auto sig = 0.05; // Reject H0 (generator is correct) with a significance level of 5%
    const auto numGraphs = 100; // Compute 100 graphs
    const auto corrected_sig = sig / n / (n - 1) * 2.0; // Bonferroni correction as we carry out (n over 2) independent trails

    auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
    std::vector<double> radii, angles;
    std::tie(radii, angles) = hypergirgs::sampleRadiiAndAngles(n, alpha, R, sampleSeed + n);

    const auto counts = count_edges(radii, angles, T, R, edgesSeed + n, numGraphs);
    const auto violations = count_violations(radii, angles, T, R, counts, numGraphs, corrected_sig);
    ASSERT_EQ(violations, 0);
}

static std::vector< std::tuple<unsigned, unsigned, double, double> > params({
    {1000, 10, 0.75, 0.5}, {900, 100, 0.75, 0.5}, {800, 10, 0.6, 0.5}, {700, 10, 0.75, 0.9}
});
INSTANTIATE_TEST_SUITE_P(Params, EdgeProbabilities_test,
                         ::testing::ValuesIn(params.begin(), params.end()));