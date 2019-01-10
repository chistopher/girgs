
#pragma once

#include <vector>
#include <cmath>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {

HYPERGIRGS_API double calculateRadius(int n, double alpha, double T, int deg);
HYPERGIRGS_API double hyperbolicDistance(double r1, double phi1, double r2, double phi2);

HYPERGIRGS_API std::vector<double> sampleRadii(int n, double alpha, double R, int seed);
HYPERGIRGS_API std::vector<double> sampleAngles(int n, int seed);

} // namespace hypergirgs
