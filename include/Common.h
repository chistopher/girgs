
#pragma once

#include <vector>


struct Node {
    std::vector<double> coord;
    double weight;
    int index; // for debug // TODO maybe remove?
    std::vector<Node*> edges;
};