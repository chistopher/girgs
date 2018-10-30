
#pragma once

#include <vector>

#include <Common.h>


class Generator
{
public:
   void generateGIRG(unsigned int dimension, const std::vector<double> &weights, double alpha, double c, int seed);

protected:


    std::vector<Node> m_graph;
};
