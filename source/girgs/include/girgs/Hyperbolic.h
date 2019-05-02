
#include <cmath>

#define M_PI       3.14159265358979323846   // pi

namespace girgs {

static double radiusToGirgWeight(double r, double R) { return std::exp((R - r) / 2); }
static double girgWeightToRadius(double w, double R, double scaling = 1.0) { return R - 2 * std::log(w / scaling); }

static double angleToGirgPosition(double angle) { return angle / 2 / M_PI; }
static double girgPositionToAngle(double position) { return position * 2 * M_PI; }

} // namespace girgs
