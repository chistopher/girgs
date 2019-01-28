#pragma once

#include <vector>
#include <utility>
#include <random>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>
#include <hypergirgs/Point.h>


namespace hypergirgs {


template <typename EdgeCallback>
class HyperbolicTree
{
public:

    HyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false);

    void generate(int seed);

protected:


    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);
    void visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
                                       unsigned int first_parallel_level, std::vector<std::vector<unsigned int>>& parallel_calls);

    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    unsigned int partitioningBaseLevel(double r1, double r2); // takes lower bound on radius for two layers

    double connectionProb(double dist); // connection probability with respect to hyperbolic distance


protected:
    EdgeCallback& m_edgeCallback;
    const bool m_profile;

    const size_t m_n; ///< number of nodes

    const double m_coshR; ///< = cosh(R)

    const double m_T; ///< temperature
    const double m_R; ///< radius

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::vector<RadiusLayer> m_radius_layers; ///< data structure holding the points

    std::vector<std::vector<std::pair<unsigned int, unsigned int> > > m_layer_pairs;

    std::vector<hypergirgs::default_random_engine> m_gens; ///< random generators for each thread
    std::vector<std::uniform_real_distribution<>> m_dists; ///< random distributions for each thread

#ifndef NDEBUG
    long long m_type1_checks; ///< number of node pairs per thread that are checked via a type 1 check
    long long m_type2_checks; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};

template <typename EdgeCallback>
inline HyperbolicTree<EdgeCallback> makeHyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false) {
    return {radii, angles, T, R, edgeCallback, profile};
}

} // namespace hypergirgs

#include <hypergirgs/HyperbolicTree.inl>
