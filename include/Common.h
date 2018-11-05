
#pragma once

#include <vector>

struct Node {
    std::vector<double> coord;
    double weight;
    int index;
    std::vector<Node*> edges;
};


std::vector<double> generateWeights(unsigned int n, double beta, int seed);
