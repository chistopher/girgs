
#pragma once

#include <vector>

struct Node {
    std::vector<double> coord;
    double weight;
    int index;
    std::vector<Node*> edges;
};


std::vector<double> generateWeights(unsigned int n, double beta, int seed);

// max over the torus distance in all dimensions
double distance(const std::vector<double>& a, const std::vector<double>& b);