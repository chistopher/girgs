
#pragma once

#include <vector>
#include <random>
#include <utility>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {

using default_random_engine = std::mt19937_64;

HYPERGIRGS_API double calculateRadius(int n, double alpha, double T, double deg);
HYPERGIRGS_API double calculateRadiusLikeNetworKit(int n, double alpha, double T, double deg);

HYPERGIRGS_API std::vector<double> sampleRadii(int n, double alpha, double R, int seed, bool parallel = true);
HYPERGIRGS_API std::vector<double> sampleAngles(int n, int seed, bool parallel = true);

/// If both, radii and angles, are to be sampled prefer this function of sampleRadii() and sampleAngles() for performance and quality reasons.
HYPERGIRGS_API std::pair<std::vector<double>, std::vector<double> > sampleRadiiAndAngles(int n, double alpha, double R, int seed, bool parallel = true);


HYPERGIRGS_API std::vector<std::pair<int, int> > generateEdges(std::vector<double>& radii, std::vector<double>& angles, double T, double R, int seed = 0);

} // namespace hypergirgs
