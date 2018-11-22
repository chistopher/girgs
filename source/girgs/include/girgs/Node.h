
#pragma once

#include <vector>

#include <girgs/girgs_api.h>


struct GIRGS_API Node {
    std::vector<double> coord;
    double weight;
    int index;
    std::vector<Node*> edges;
};


// max over the torus distance in all dimensions TODO find place for this function
GIRGS_API double distance(const std::vector<double>& a, const std::vector<double>& b);
