
#pragma once

#include <vector>

struct Node {
    std::vector<double> coord;
    double weight;
    int index;
    std::vector<Node*> edges;
};


// max over the torus distance in all dimensions TODO find place for this function
double distance(const std::vector<double>& a, const std::vector<double>& b);