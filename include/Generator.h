
#pragma once

#include <vector>
#include <string>

#include <Common.h>


// a generator that samples a GIRG in expected linear time and provides access to the result
// can be used multiple times
// if multiple graphs are generated with the same generator instance, the generation process implies the deletion of the previous graph
// accessors to graph data always refer to the last generated graph
// weightSeed and position seed should not be equal!
class Generator
{
public:

    // init parameter
    void setWeights(const std::vector<double>& weights);
    void setWeights(int n, double ple, int weightSeed);
    void setPositions(const std::vector<std::vector<double>>& positions);
    void setPositions(int n, int dimension, int positionSeed);

    double scaleWeights(int desiredAvgDegree, int dimension, double alpha);

    // generate
    void generate(double alpha, int samplingSeed);
    void generateTreshold();
    // convenience
    const std::vector<Node>& generate(int n, int dimension, double ple, double alpha, int desiredAvgDegree, int weightSeed, int positionSeed, int samplingSeed);


    // access results
    const std::vector<Node>& graph() const {return m_graph;}
    double avg_degree() const;
    void saveDot(std::string file) const;
    // copy internal data for access
    std::vector<double> weights() const;
    std::vector<std::vector<double>> positions() const;

    // helper
    // currently works only for alpha = std::numeric_limits<double>::infinity()
    double estimateWeightScaling(const std::vector<double>& weights, int desiredAvgDegree, int dimension, double alpha) const;

protected:

    std::vector<Node> m_graph;
};
