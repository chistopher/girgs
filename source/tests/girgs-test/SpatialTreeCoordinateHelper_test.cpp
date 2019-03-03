
#include <random>

#include <gmock/gmock.h>

#include <girgs/SpatialTreeCoordinateHelper.h>


using namespace std;
using namespace girgs;


class SpatialTreeCoordinateHelper_test: public testing::Test
{
protected:
    SpatialTreeCoordinateHelper_test() : b1(12), b2(6), b3(4), b4(3) {}

    SpatialTreeCoordinateHelper<1> b1;
    SpatialTreeCoordinateHelper<2> b2;
    SpatialTreeCoordinateHelper<3> b3;
    SpatialTreeCoordinateHelper<4> b4;
};


template<unsigned int D>
void testTreeStructure(SpatialTreeCoordinateHelper<D>& tree) {

    // test consistency of cell index up to level 12 / D
    const auto max_level = tree.levels();

    // check the number of children on all layers and helper functions
    auto cells = tree.firstCellOfLevel(max_level+1);
    auto sumCells = 1; // root in level 0
    auto numChildren = vector<int>(cells, 0);
    for(auto l=1u; l<=max_level; ++l) {
        EXPECT_EQ(sumCells, tree.firstCellOfLevel(l));
        EXPECT_EQ(sumCells, tree.parent(tree.firstCellOfLevel(l+1)));
        EXPECT_EQ(sumCells, tree.firstChild(tree.firstCellOfLevel(l-1)));
        EXPECT_EQ(sumCells, tree.firstCellOfLevel(l+1) - tree.numCellsInLevel(l));
        sumCells += tree.numCellsInLevel(l);
        for(auto i = tree.firstCellOfLevel(l); i<tree.firstCellOfLevel(l+1); ++i) {
            numChildren.at(tree.parent(i)) += 1; // use at to get bounds check
        }
    }
    EXPECT_EQ(cells, sumCells);

    // check that all but the last level have correct number of children
    for(auto i=0; i<tree.firstCellOfLevel(max_level); ++i)
        EXPECT_EQ(numChildren[i], tree.numChildren());
    for(auto i=tree.firstCellOfLevel(max_level); i<cells; ++i)
        EXPECT_EQ(numChildren[i], 0);
}

template<unsigned int D>
void testCoordMapping(SpatialTreeCoordinateHelper<D>& helper) {

    // generate some points and check their cells on all levels
    mt19937 gen(1337);
    uniform_real_distribution<> dist(0.0, 1.0);
    for(auto i=0; i<100; ++i) {

        // generate random point in [0..1)^D
        auto point = std::array<double, D>();
        for(auto d=0u; d<D; ++d)
            point[d] = dist(gen);

        // compute containing cell in all levels and check if point is in their bounds
        auto containingCells = vector<unsigned int>(helper.levels());
        containingCells[0] = 0;
        for(auto l=1u; l<helper.levels(); ++l){
            containingCells[l] = helper.cellForPoint(point, l) + helper.firstCellOfLevel(l);
            auto cell = containingCells[l];
            auto bounds = helper.bounds(cell, l);
            for(auto d=0u; d<D; ++d){
                // point is in cell bounds
                EXPECT_LE(bounds[d].first, point[d]);
                EXPECT_LT(point[d], bounds[d].second);
            }
        }

        // check that all containing cells have the same parents
        for(auto l=1; l<helper.levels(); ++l)
            EXPECT_EQ(containingCells[l-1], helper.parent(containingCells[l]));
    }
}


TEST_F(SpatialTreeCoordinateHelper_test, testTreeStructure)
{
    testTreeStructure(b1);
    testTreeStructure(b2);
    testTreeStructure(b3);
    testTreeStructure(b4);
}


TEST_F(SpatialTreeCoordinateHelper_test, testCoordMapping)
{
    testCoordMapping(b1);
    testCoordMapping(b2);
    testCoordMapping(b3);
    testCoordMapping(b4);
}


TEST_F(SpatialTreeCoordinateHelper_test, testTouching)
{
    // TODO write better test for touching
    auto a = b1.firstCellOfLevel(5);
    auto b = b1.firstCellOfLevel(5) + b1.numCellsInLevel(5) - 1;
    EXPECT_TRUE(b1.touching(a,b,5));
}


TEST_F(SpatialTreeCoordinateHelper_test, testDistance)
{
    // TODO write better tests for dist of cells
    EXPECT_EQ(b1.dist(10, 11, 3), 0);
    EXPECT_EQ(b1.dist(10, 12, 3), (1.0/8)* 1);
    EXPECT_EQ(b1.dist(10, 13, 3), (1.0/8)* 2);
    EXPECT_EQ(b1.dist(10, 14, 3), (1.0/8)* 3);

    EXPECT_EQ(b2.dist(10,  8, 2), (1.0/4)* 1);
    EXPECT_EQ(b2.dist(10,  9, 2), 0);
    EXPECT_EQ(b2.dist(10, 14, 2), (1.0/4)* 1);
    EXPECT_EQ(b2.dist(10, 17, 2), (1.0/4)* 1);
}
