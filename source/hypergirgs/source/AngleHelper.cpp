#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/AngleHelper.h>

#include <cassert>
#include <cmath>
#include <algorithm>

namespace hypergirgs {


std::pair<double, double> AngleHelper::bounds(unsigned int cell, unsigned int level) {
    auto diameter = 2*PI / (1<<level);
    auto localIndex = cell - firstCellOfLevel(level);
    return {localIndex*diameter, (localIndex+1) * diameter};
}

bool AngleHelper::touching(unsigned int cellA, unsigned int cellB, unsigned int level) {
    auto mm = std::minmax(cellA,cellB);
    auto diff = mm.second - mm.first;
    return diff<=1 || diff == numCellsInLevel(level) - 1;
}

double AngleHelper::dist(unsigned int cellA, unsigned int cellB, unsigned int level) {
    auto mm = std::minmax(cellA,cellB);
    auto diff = mm.second - mm.first;
    return (diff <= 1) ? 0 : (diff-1) * 2.0*PI / (1<<level);
}


} // namespace hypergirgs
