
#pragma once

#include <vector>
#include <utility>
#include <random>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>
#include <hypergirgs/Point.h>

#include <hypergirgs/hypergirgs_api.h>


namespace hypergirgs {


class HYPERGIRGS_API HyperbolicTree
{
public:

    HyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R);

    std::vector<std::pair<int,int>> generate(int seed);

protected:


    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);

    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    unsigned int partitioningBaseLevel(double r1, double r2); // takes lower bound on radius for two layers


protected:

    std::vector<std::pair<int,int>> m_result;

    const size_t m_n; ///< number of nodes

    const double m_coshR; ///< = cosh(R)

    const double m_T;
    const double m_R;

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::vector<RadiusLayer> m_radius_layers;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs;

    hypergirgs::default_random_engine m_gen; ///< random generator
    std::uniform_real_distribution<> m_dist; ///< random distribution

#ifndef NDEBUG
    long long m_type1_checks; ///< number of node pairs per thread that are checked via a type 1 check
    long long m_type2_checks; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};


} // namespace hypergirgs
