
#include <iostream>
#include <random>

#include <SpatialTree.h>


using namespace std;

void test(bool cond){
    if(!cond){
        cout << "assertion failed" << endl;
        exit(1);
    }
}


template<class T>
void testTreeStructure(T& tree) {

    // check the number of children on all layers and helper functions
    auto cells = T::numCells;
    auto sumCells = 1; // root in level 0
    auto numChildren = array<int, T::numCells>();
    numChildren.fill(0);
    for(auto l=1u; l<T::level; ++l) {
        test(sumCells == T::firstCellOfLevel(l));
        test(sumCells == T::parent(T::firstCellOfLevel(l+1)));
        test(sumCells == T::firstChild(T::firstCellOfLevel(l-1)));
        test(sumCells + T::numCellsInLevel(l) == T::firstCellOfLevel(l+1));
        sumCells += T::numCellsInLevel(l);
        for(auto i = T::firstCellOfLevel(l); i<T::firstCellOfLevel(l+1); ++i) {
            numChildren.at(T::parent(i)) += 1; // use array::at to get bounds check
        }
    }
    test(cells == sumCells);

    // check that all but the last level have correct number of children
    for(auto i=0; i<T::firstCellOfLevel(T::level-1); ++i)
        test(numChildren[i] == T::numChildren);
    for(auto i=T::firstCellOfLevel(T::level-1); i<cells; ++i)
        test(numChildren[i] == 0);
}

template<class T>
void testCoordMapping(T& tree) {

    // generate some points and check their cells on all levels
    mt19937 gen(1337);
    uniform_real_distribution<> dist(0.0, 1.0);
    for(auto i=0; i<100; ++i) {

        // generate random point in [0..1)^D
        array<double, T::dimension> point;
        for(auto d=0; d<T::dimension; ++d)
            point[d] = dist(gen);

        // compute containing cell in all levels and check if point is in their bounds
        array<unsigned int, T::level> containingCells;
        containingCells[0] = 0;
        for(auto l=1u; l<T::level; ++l){
            containingCells[l] = tree.cellForPoint(point, l);
            auto cell = containingCells[l];
            auto bounds = tree.bounds(cell, l);
            for(auto d=1u; d<T::dimension; ++d){
                test(bounds[d].first <= point[d] && point[d] < bounds[d].second);
            }
        }

        // check that all containing cells have the same parents
        for(auto l=1; l<T::level; ++l)
            test(containingCells[l-1] == T::parent(containingCells[l]));
    }
}


int main(int argc, char* argv[]) {

    constexpr auto a1 = SpatialTree<1,7>();
    constexpr auto a2 = SpatialTree<2,6>();
    constexpr auto a3 = SpatialTree<3,5>();
    constexpr auto a4 = SpatialTree<4,4>();


    testTreeStructure(a1);
    testTreeStructure(a2);
    testTreeStructure(a3);
    testTreeStructure(a4);

    testCoordMapping(a1);
    testCoordMapping(a2);
    testCoordMapping(a3);
    testCoordMapping(a4);

    cout << "all tests passed." << endl;
    return 0;
}
