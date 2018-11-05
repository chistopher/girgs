
#include <Generator.h>

#include <fstream>
#include <iomanip>

#include <SpatialTree.h>


void Generator::generateGIRG(unsigned int dimension, const std::vector<double> &weights, double alpha, double c, int seed) {

    if(0 == dimension || dimension > 5) throw std::logic_error("dimension to high");
    m_graph =
            (dimension == 1) ? SpatialTree<1>().generateGraph(weights, alpha, c, seed) :
            (dimension == 2) ? SpatialTree<2>().generateGraph(weights, alpha, c, seed) :
            (dimension == 3) ? SpatialTree<3>().generateGraph(weights, alpha, c, seed) :
            (dimension == 4) ? SpatialTree<4>().generateGraph(weights, alpha, c, seed) :
            (dimension == 5) ? SpatialTree<5>().generateGraph(weights, alpha, c, seed) :
            std::vector<Node>();
}

double Generator::avg_degree() const {
    auto edges = 0.0;
    for(auto& each : graph())
        edges += each.edges.size();
    return edges / m_graph.size();
}


void Generator::saveDot(std::string file) const {
    if(m_graph.empty()){
        std::cout << "no graph generated" << std::endl;
        return;
    }

    auto f = std::ofstream(file);
    f << "graph girg {\n\toverlap=scale;\n\n";
    for(auto& each : m_graph){
        f   << '\t' << each.index << " [label=\""
            << std::setprecision(2) << std::fixed << each.weight << std::defaultfloat << std::setprecision(6)
            << "\", pos=\"";
        for(auto d=0u; d<each.coord.size(); ++d)
            f << (d==0 ? "" : ",") << each.coord[d];
        f << "\"];\n";
    }
    f << '\n';
    for(auto& each : m_graph){
        for(auto neighbor : each.edges)
            if(each.index < neighbor->index)
                f << '\t' << each.index << "\t-- " << neighbor->index << ";\n";
    }
    f << "}\n";
}

