
#pragma once

#include <vector>
#include <utility>
#include <random>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>

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
    const std::vector<double> m_radii;
    const std::vector<double> m_angles;

    double m_pre_R;
    std::vector<std::pair<std::pair<double, double>,double>> m_pre_coord; ///< pre-computations for distance

    const double m_T;
    const double m_R;

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::vector<RadiusLayer> m_radius_layers;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs;
   
    std::mt19937 m_gen; ///< random generator
    std::uniform_real_distribution<> m_dist; ///< random distribution

#ifndef NDEBUG
    long long m_type1_checks; ///< number of node pairs per thread that are checked via a type 1 check
    long long m_type2_checks; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};


} // namespace hypergirgs
