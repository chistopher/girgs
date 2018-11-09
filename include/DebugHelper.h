
#pragma once

#include <iostream>
#include <vector>
#include <random>


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

std::vector<double> generatePowerLawWeights(unsigned int n, double lower, double upper, double beta, int seed) {
    auto gen = std::mt19937(seed >= 0 ?  seed : std::random_device()());
    std::uniform_real_distribution<> dist(0,1);
    auto weights = std::vector<double>(n);
    for(auto i=0; i<n; ++i)
        weights[i] = std::pow(
                (std::pow(upper,beta+1)-std::pow(lower,beta+1))*dist(gen) + std::pow(lower,beta+1),
                1.0/(beta+1) );
    return weights;
}