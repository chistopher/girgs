
#include <random>

#include <gmock/gmock.h>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/Generator.h>
#include <hypergirgs/Point.h>


using namespace std;
using namespace hypergirgs;


class AngleHelper_test: public testing::Test
{
protected:
    AngleHelper_test()
    : n(100)
    , angles(sampleAngles(n, 12))
    , radii(sampleRadii(n, 0.75, calculateRadius(n, 0.75, 0, 10), 15))
    {
    }

    unsigned int n;
    vector<double> angles;
    vector<double> radii;
};


TEST_F(AngleHelper_test, testCellInsertion)
{
    // test that all points are inserted in the correct cells
    const auto maxLevel = 10u;
    for (unsigned int level = 0; level < maxLevel; ++level) {
        for(auto angle : angles) {
            auto cell = AngleHelper::cellForPoint(angle, level) + AngleHelper::firstCellOfLevel(level);
            auto bounds = AngleHelper::bounds(cell, level);
            ASSERT_GE(angle, bounds.first);
            ASSERT_LT(angle, bounds.second);
        }
    }
}


TEST_F(AngleHelper_test, testTouching)
{
    auto level = 5u;

    auto a1 = 2.0;
    auto a2 = 4.0;

    auto cell1 = AngleHelper::cellForPoint(a1, level) + AngleHelper::firstCellOfLevel(level);
    auto cell2 = AngleHelper::cellForPoint(a2, level) + AngleHelper::firstCellOfLevel(level);
    ASSERT_FALSE(AngleHelper::touching(cell1, cell2, level));
}


TEST_F(AngleHelper_test, testDistance)
{
    auto level = 7u;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            auto cell1 = AngleHelper::cellForPoint(angles[i], level) + AngleHelper::firstCellOfLevel(level);
            auto cell2 = AngleHelper::cellForPoint(angles[j], level) + AngleHelper::firstCellOfLevel(level);
            auto actualDist = hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]);
            auto cellDist = AngleHelper::dist(cell1, cell2, level);
            ASSERT_LE(cellDist, actualDist);
        }
    }
}
