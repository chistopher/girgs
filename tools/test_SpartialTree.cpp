
#include <iostream>
#include <random>

#include <DebugHelper.h>
#include <SpatialTree.h>
#include <SpatialTreeCoordinateHelper.h>


using namespace std;

void test(bool cond){
    if(!cond){
        cout << "assertion failed" << endl;
        exit(1);
    }
}


template<unsigned int D>
void testTreeStructure(SpatialTree<D>& tree) {

    // test consistency of cell index up to level 12 / D
    const auto max_level = 12 / D;

    // check the number of children on all layers and helper functions
    auto cells = tree.firstCellOfLevel(max_level+1);
    auto sumCells = 1; // root in level 0
    auto numChildren = vector<int>(cells, 0);
    for(auto l=1u; l<=max_level; ++l) {
        test(sumCells == tree.firstCellOfLevel(l));
        test(sumCells == tree.parent(tree.firstCellOfLevel(l+1)));
        test(sumCells == tree.firstChild(tree.firstCellOfLevel(l-1)));
        test(sumCells + tree.numCellsInLevel(l) == tree.firstCellOfLevel(l+1));
        sumCells += tree.numCellsInLevel(l);
        for(auto i = tree.firstCellOfLevel(l); i<tree.firstCellOfLevel(l+1); ++i) {
            numChildren.at(tree.parent(i)) += 1; // use at to get bounds check
        }
    }
    test(cells == sumCells);

    // check that all but the last level have correct number of children
    for(auto i=0; i<tree.firstCellOfLevel(max_level); ++i)
        test(numChildren[i] == tree.numChildren);
    for(auto i=tree.firstCellOfLevel(max_level); i<cells; ++i)
        test(numChildren[i] == 0);
}

template<unsigned int D>
void testCoordMapping(SpatialTreeCoordinateHelper<D>& helper) {

    // generate some points and check their cells on all levels
    mt19937 gen(1337);
    uniform_real_distribution<> dist(0.0, 1.0);
    for(auto i=0; i<100; ++i) {

        // generate random point in [0..1)^D
        auto point = vector<double>(D, 0.0);
        for(auto d=0u; d<D; ++d)
            point[d] = dist(gen);

        // compute containing cell in all levels and check if point is in their bounds
        auto containingCells = vector<unsigned int>(helper.levels());
        containingCells[0] = 0;
        for(auto l=1u; l<helper.levels(); ++l){
            containingCells[l] = helper.cellForPoint(point, l);
            auto cell = containingCells[l];
            auto bounds = helper.bounds(cell, l);
            for(auto d=0u; d<D; ++d){
                test(bounds[d].first <= point[d] && point[d] < bounds[d].second);
            }
        }

        // check that all containing cells have the same parents
        for(auto l=1; l<helper.levels(); ++l)
            test(containingCells[l-1] == helper.parent(containingCells[l]));
    }
}



int main(int argc, char* argv[]) {

    auto a1 = SpatialTree<1>();
    auto a2 = SpatialTree<2>();
    auto a3 = SpatialTree<3>();
    auto a4 = SpatialTree<4>();

    testTreeStructure(a1);
    testTreeStructure(a2);
    testTreeStructure(a3);
    testTreeStructure(a4);

    auto b1 = SpatialTreeCoordinateHelper<1>(12);
    auto b2 = SpatialTreeCoordinateHelper<2>(6);
    auto b3 = SpatialTreeCoordinateHelper<3>(4);
    auto b4 = SpatialTreeCoordinateHelper<4>(3);

    testCoordMapping(b1);
    testCoordMapping(b2);
    testCoordMapping(b3);
    testCoordMapping(b4);

    // TODO write better test for touching
    auto a = a1.firstCellOfLevel(5);
    auto b = a1.firstCellOfLevel(5) + a1.numCellsInLevel(5) - 1;
    test(b1.touching(a,b,5));

    // TODO write better tests for dist of cells
    test(b1.dist(10, 11, 3) == 0);
    test(b1.dist(10, 12, 3) == (1.0/8)* 1);
    test(b1.dist(10, 13, 3) == (1.0/8)* 2);
    test(b1.dist(10, 14, 3) == (1.0/8)* 3);

    test(b2.dist(10, 8, 2) == (1.0/4)* 1);
    test(b2.dist(10, 9, 2) == 0);
    test(b2.dist(10, 14, 2) == (1.0/4)* 1);
    test(b2.dist(10, 17, 2) == (1.0/4)* 1);


    cout << "all tests passed." << endl;
    return 0;
}
