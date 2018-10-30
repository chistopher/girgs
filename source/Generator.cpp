
#include <Generator.h>

#include <SpatialTree.h>


void
Generator::generateGIRG(unsigned int dimension, const std::vector<double> &weights, double alpha, double c, int seed) {

    if(0 == dimension || dimension > 5) throw std::logic_error("dimension to high");
    m_graph =
            (dimension == 1) ? SpatialTree<1>().generateGraph(weights, alpha, c, seed) :
            (dimension == 2) ? SpatialTree<2>().generateGraph(weights, alpha, c, seed) :
            (dimension == 3) ? SpatialTree<3>().generateGraph(weights, alpha, c, seed) :
            (dimension == 4) ? SpatialTree<4>().generateGraph(weights, alpha, c, seed) :
            (dimension == 5) ? SpatialTree<5>().generateGraph(weights, alpha, c, seed) :
            std::vector<Node>();
}
