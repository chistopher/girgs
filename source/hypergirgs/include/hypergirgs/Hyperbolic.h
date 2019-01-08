
#include <vector>
#include <string>
#include <cmath>

#include <hypergirgs/hypergirgs_api.h>

namespace hypergirgs {

HYPERGIRGS_API double calculateRadius(int n, double alpha, double T, int deg);
HYPERGIRGS_API double hyperbolicDistance(double r1, double phi1, double r2, double phi2);

HYPERGIRGS_API std::vector<double> sampleRadii(int n, double alpha, double R, int weightSeed);
HYPERGIRGS_API std::vector<double> sampleAngles(int n, int positionSeed);

double radiusToGirgWeight(double r, double R) { return std::exp((R - r) / 2); }
double girgWeightToRadius(double w, double R, double scaling = 1.0) { return R - 2 * std::log(w / scaling); }

double angleToGirgPosition(double angle) { return angle / 2 / M_PI; }
double girgPositionToAngle(double position) { return position * 2 * M_PI; }

} // namespace hypergirgs
