
#pragma once

#include <vector>
#include <string>

#include <Common.h>


// a generator that samples a GIRG in expected linear time and provides access to the result
// can be used multiple times
// if multiple graphs are generated with the same generator instance, the generation process implies the deletion of the previous graph
// accessors to graph data always refer to the last generated graph
class Generator
{
public:
   void generateGIRG(unsigned int dimension, const std::vector<double> &weights, double alpha, double c, int seed);

   const std::vector<Node>& graph() const {return m_graph;}

   double avg_degree() const;

   void saveDot(std::string file) const;

protected:


    std::vector<Node> m_graph;
};
