
#include <vector>
#include <string>
#include <cmath>

#include <girgs/Node.h>
#include <girgs/girgs_api.h>

namespace girgs {

GIRGS_API double calculateRadius(int n, double alpha, double T, int deg);
GIRGS_API double hyperbolicDistance(double r1, double phi1, double r2, double phi2);

GIRGS_API std::vector<double> sampleRadii(int n, double alpha, double R, int weightSeed);
GIRGS_API std::vector<double> sampleAngles(int n, int positionSeed);

GIRGS_API double radiusToGirgWeight(double r, double R) { return std::exp((R - r) / 2); }
GIRGS_API double girgWeightToRadius(double w, double R, double scaling = 1.0) { return R - 2 * std::log(w / scaling); }

GIRGS_API double angleToGirgPosition(double angle) { return angle / 2 / M_PI; }
GIRGS_API double girgPositionToAngle(double position) { return position * 2 * M_PI; }

} // namespace girgs
