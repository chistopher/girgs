
#pragma once

#include <iostream>

template<class T>
void printInfo(const T & tree) {
    using namespace std;
    cout << "dimension: " << tree.dimension << endl;
    cout << "level: " << tree.level << endl;
    cout << "cells: " << tree.numCells << endl;
    cout << "children per cell: " << tree.numChildren << endl;
    cout << endl;
    cout << "index ranges of levels" << endl;
    cout << "l\tcells\tindexRange" << endl;
    for(int i=0; i<tree.level; ++i)
      cout << i << '\t' << tree.numCellsInLevel(i)
          << "\t[" << tree.firstCellOfLevel(i) << ".." << tree.firstCellOfLevel(i) + tree.numCellsInLevel(i) - 1<< "]\n";

    cout << endl;
    cout << "cell / parent / first child of parent / i-th child of parent" << endl;
    for(auto l=1; l<tree.level; ++l) {
        cout << endl << "Level: " << l << endl;
        for(int i=0; i<tree.numCellsInLevel(l); ++i){
            auto child = i + tree.firstCellOfLevel(l);
            auto parent = tree.parent(child);
            cout << child << "\t" << parent << "\t" << tree.firstChild(parent) << "\t" << ((child-1)&(tree.numChildren-1)) << endl;
        }
    }


    return;
}
