
#include <cassert>
#include <cmath>
#include <algorithm>

#include <hypergirgs/AngleHelper.h>



namespace hypergirgs {


std::pair<double, double> AngleHelper::bounds(unsigned int cell, unsigned int level) {
    auto diameter = 2*M_PI / (1<<level);
    auto localIndex = cell - firstCellOfLevel(level);
    return {localIndex*diameter, (localIndex+1) * diameter};
}

unsigned int AngleHelper::cellForPoint(double angle, unsigned int targetLevel) {
    return static_cast<unsigned int>(angle/2/M_PI * numCellsInLevel(targetLevel));
}

bool AngleHelper::touching(unsigned int cellA, unsigned int cellB, unsigned int level) {
    auto mm = std::minmax(cellA,cellB);
    auto diff = mm.second - mm.first;
    return diff<=1 || diff == numCellsInLevel(level) - 1;
}

double AngleHelper::dist(unsigned int cellA, unsigned int cellB, unsigned int level) {
    auto mm = std::minmax(cellA,cellB);
    auto diff = mm.second - mm.first;
    return (diff <= 1) ? 0 : (diff-1) * 2*M_PI / (1<<level);
}


} // namespace hypergirgs
