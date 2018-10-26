
#pragma once

#include <array>
#include <vector>


template<unsigned int D>
struct Node {
    std::array<double, D> coord;
    double weight;
    int index; // for debug // TODO maybe remove?
    std::vector<Node*> edges;
};