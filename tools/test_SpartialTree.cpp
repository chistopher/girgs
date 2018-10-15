
#include <iostream>
#include <SpatialTree.h>


using namespace std;

int main(int argc, char* argv[])
{
    //auto d1l1_points =
    auto d1l4 = SpatialTree<1,4>();

    auto a = SpatialTree<3,4>();

    //cout << a.numCells << endl;

    int sumCells = 0;
    for(int i=0; i<4; ++i)
        sumCells += a.numCellsInLevel(i);
    //cout << sumCells << endl;

    cout << "dimension: " << a.dimension << endl;
    cout << "level: " << a.level << endl;
    cout << "cells: " << a.numCells << endl;
    cout << "children per cell: " << a.numChildren << endl;
    cout << endl;
    cout << "index ranges of levels" << endl;
    cout << "l\tcells\tindexRange" << endl;
    for(int i=0; i<a.level; ++i)
        cout << i << '\t' << a.numCellsInLevel(i)
            << "\t[" << a.firstCellOfLevel(i) << ".." << a.firstCellOfLevel(i) + a.numCellsInLevel(i) - 1<< "]\n";

    cout << endl;
    cout << "info for all cells of level 2" << endl;
    cout << "cell / parent / first child of parent / i-th child of parent" << endl;
    for(int i=0; i<a.numCellsInLevel(2); ++i){
        auto child = i + a.firstCellOfLevel(2);
        auto parent = a.parent(child);
        cout << child << "\t" << parent << "\t" << a.firstChild(parent) << "\t" << ((child-1)&(a.numChildren-1)) << endl;
    }

    // cout << a.pointsInCell()

    return 0;
}
