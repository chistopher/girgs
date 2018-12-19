
#pragma once

#include <vector>
#include <string>
#include <functional>

#include <girgs/girgs_api.h>
#include <girgs/Node.h>


namespace girgs {


/**
 * @brief
 *  A generator that samples a GIRG in expected linear time and provides access to the result.
 *  The generator can be used for multiple GIRGs.
 *  If multiple graphs are generated with the same generator instance, the generation process implies the deletion of the previous graph.
 *  Accessors to graph data always refer to the last generated graph.
 *  WeightSeed and position seed should not be equal!
 */
class GIRGS_API Generator
{
public:

    /**
     * @brief
     *  Set new weights explicitly.
     *
     * @param weights
     *  Weights for the generation process. Size of the vector should match with positions.
     *  The behaviour is undefined if the weights are negative or do not follow a power law distribution.
     */
    void setWeights(const std::vector<double>& weights);

    /**
     * @brief
     *  Set new weights implicitly. The weights are sampled according to a power law distribution between [1, n)
     *
     * @param n
     *  The size of the graph. Should match with size of positions.
     * @param ple
     *  The power law exponent to sample the new weights. Should be negative and approximately -2.0 to -3.0.
     * @param weightSeed
     *  A seed for weight sampling. Should not be equal to the position seed.
     */
    void setWeights(int n, double ple, int weightSeed);

    /**
     * @brief
     *  Set new positions explicitly.
     *
     * @param positions
     *  Explicit positions for the underlying geometry. The coordinates for all nodes should have the same size.
     */
    void setPositions(const std::vector<std::vector<double>>& positions);

    /**
     * @brief
     *  Samples d dimensional coordinates for n points on a torus \f$[0,1)^d\f$.
     * @param n
     *  Size of the graph.
     * @param dimension
     *  Dimension of the geometry.
     * @param positionSeed
     *  Seed to sample the positions.
     */
    void setPositions(int n, int dimension, int positionSeed);

    /**
     * @brief
     *  Scales all weights so that the expected average degree equals desiredAvgDegree.
     *  Implemented as binary search over an estimation function.
     *
     * @bug
     *  For \f$\alpha > 10\f$ we use the estimation for threshold graphs due to numerical difficulties.
     *  This leads to slightly higher degrees than desired.
     *  Also I experienced inaccurate results for \f$9 \leq \alpha < 10\f$.
     *
     * @param desiredAvgDegree
     *  The desired average degree
     * @param dimension
     *  Dimension of the underlying geometry. Should equal the dimensionality of the positions.
     * @param alpha
     *  Parameter of the algorithm. Should be the same as for the generation process.
     *
     * @return
     *  The scaling s applied to all weights.
     *  The constant c hidden in the theta of the edge probabilities is \f$s^\alpha\f$ for \f$\alpha < \infty\f$
     *  and \f$s^{1/d}\f$ in the threshold case.
     */
    double scaleWeights(double desiredAvgDegree, int dimension, double alpha);

    /**
     * @brief
     *  Samples edges according to the current weights and positions.
     *  An edge between node u and v is formed with probability \f$ \left(\frac{w_u w_v / W}{|| x_u - x_v ||^d}\right)^\alpha \f$ or 1.0 if the term exceeds 1.0.
     *
     * @pre Weights and positions must be set and have equal length.
     *
     * @param alpha
     *  Parameter of the algorithm. Using infinity results in a deterministic threshold graph
     *  and is equivalent to calling generateThreshold().
     * @param samplingSeed
     *  A seed used for the generation process.
     */
    void generate(double alpha, int samplingSeed);

    /**
     * @brief
     *  Convenience method that sets weights and positions, scales weights, and samples the edges.
     *  This enables code like:
     *
     *  `auto graph = girgs::Generator().generate(...);`
     *
     *  or
     *
     *  `for(auto& node : girgs::Generator().generate(...))`
     *
     *  Note that the data must be moved, because it stores the edges as pointer.
     *
     * @param n
     *  Size of the graph.
     * @param dimension
     *  Dimension of the geometry.
     * @param ple
     *  The power law exponent to sample the weights. (see setWeights(int, double, int))
     * @param alpha
     *  Edge probability parameter.
     * @param desiredAvgDegree
     *  Desired average degree. (see scaleWeights(double, int, double))
     * @param weightSeed
     *  Seed to sample the weights. Should not be equal to positionSeed. (see setWeights(int, double, int))
     * @param positionSeed
     *  Seed to sample the positions. Should not be equal to weightSeed. (see setPositions(int, int, int))
     * @param samplingSeed
     *  Seed to sample the edges. (see generate(double, int))
     * @return
     *  A graph which is move constructed from the #m_graph member of the current Generator instance.
     */
    std::vector<Node> generate(int n, int dimension, double ple, double alpha, double desiredAvgDegree, int weightSeed, int positionSeed, int samplingSeed);

    /**
    *  @brief generates a threshold GIRG. (i.e. nodes u,v are connected if \f$ || x_u - x_v || < (w_u w_v / W)^{1/\alpha} \f$)
    */
    void generateThreshold();


    /**
     * @return a reference to the last sampled graph
     */
    const std::vector<Node>& graph() const {return m_graph;}

    /**
     * @return the average degree of the current graph
     */
    double avg_degree() const;

    /**
     * @return the number of edges in the current graph
     */
    unsigned int edges() const;

    /**
     * @brief
     *  Saves the graph in dot format (graphviz).
     *  The weight is saved as a label and the coordinates as a position attribute for each Node.
     *
     * @param file
     *  The name of the output file. Should end with ".dot".
     */
    void saveDot(std::string file) const;

    /**
     * @brief
     *  Saves the graph as an edge list.
     *  The first line contains the number of nodes and the number of edges.
     *  Then for each edge there is a line with the indices (zero based) of the two endpoints.
     *  All numbers are separated by spaces.
     *
     * @param file
     *  The name of the output file.
     */
    void saveEdgeList(std::string file) const;


    /**
     * @brief
     *  Provides a copy of all weights.
     *
     * @return
     *  A copy of the weights saved in the current graph (#m_graph accessed with graph() const)
     */
    std::vector<double> weights() const;

    /**
     * @brief
     *  Provides a copy of all positions.
     *
     * @return
     *  A copy of the positions saved in the current graph (#m_graph accessed with graph() const)
     */
    std::vector<std::vector<double>> positions() const;

private:
    // helper
    double estimateWeightScalingThreshold(const std::vector<double>& weights, double desiredAvgDegree, int dimension) const;
    double estimateWeightScaling(const std::vector<double>& weights, double desiredAvgDegree, int dimension, double alpha) const;

    double exponentialSearch(std::function<double(double)> f, double desiredValue, double accuracy = 0.02, double lower = 1.0, double upper = 2.0) const;

protected:

    std::vector<Node> m_graph;  ///< stores the current graph including weights and positions
};


} // namespace girgs
