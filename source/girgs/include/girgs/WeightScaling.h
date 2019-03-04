#pragma once

#include <vector>

namespace girgs {

double estimateWeightScaling(const std::vector<double> &weights, double desiredAvgDegree, int dimension, double alpha);

double estimateWeightScalingThreshold(const std::vector<double>& weights, double desiredAvgDegree, int dimension);

} // namespace girgs
