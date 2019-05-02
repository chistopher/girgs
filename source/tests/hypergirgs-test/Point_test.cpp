
#include <gmock/gmock.h>

#include <hypergirgs/Point.h>
#include <hypergirgs/Generator.h>
#include <hypergirgs/AngleHelper.h>


using namespace hypergirgs;

class Point_test: public testing::Test
{
protected:
    Point_test()
            : n(1000)
            , alpha(0.75)
            , R(calculateRadius(n, alpha, 0, 10))
            , cosh_R(std::cosh(R))
            , angles(sampleAngles(n, 12))
            , radii(sampleRadii(n, alpha, R, 15))
            , points(n)
    {
        for (int i = 0; i < n; ++i) {
            points[i] = Point(i, radii[i], angles[i], 0);
        }
    }

    unsigned int n;
    double alpha;
    double R;
    double cosh_R;
    std::vector<double> angles;
    std::vector<double> radii;
    std::vector<Point> points;
};


TEST_F(Point_test, testDistanceBelowR)
{
    for (int i = 0; i < n; ++i) {
        ASSERT_LE(0, radii[i]);
        ASSERT_LE(0, angles[i]);
        ASSERT_LT(radii[i], R);
        ASSERT_LT(angles[i], 2*PI);
        for (int j = i+1; j < n; ++j) {
            auto dist = hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]);
            if(points[i].isDistanceBelowR(points[j], cosh_R)){
                ASSERT_LT(dist, R);
            } else {
                ASSERT_GE(dist, R);
            }
        }
    }
}


TEST_F(Point_test, testHyperbolicDistance)
{
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            auto dist1 = hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]);
            auto dist2 = points[i].hyperbolicDistance(points[j]);
            auto diff = std::fabs(dist1 - dist2);
            ASSERT_LE(diff, 0.000005);
        }
    }
}
