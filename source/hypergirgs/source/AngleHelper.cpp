
#include <cassert>
#include <cmath>

#include <hypergirgs/AngleHelper.h>



namespace hypergirgs {


std::pair<double, double> AngleHelper::bounds(unsigned int cell, unsigned int level) const {

    // assert cell in correct level
    assert(firstCellOfLevel(level) <= cell && cell < firstCellOfLevel(level+1));

    auto diameter = 2*M_PI / (1<<level);
    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0u; d<D; ++d)
        result[d]= { m_coords[cell][d]*diameter, (m_coords[cell][d]+1)*diameter };
    return result;
}

unsigned int AngleHelper::cellForPoint(double angle, unsigned int targetLevel) const {
    return 0;
}

bool AngleHelper::touching(unsigned int cellA, unsigned int cellB, unsigned int level) const {
    return false;
}

double AngleHelper::dist(unsigned int cellA, unsigned int cellB, unsigned int level) const {
    return 0;
}


} // namespace hypergirgs
