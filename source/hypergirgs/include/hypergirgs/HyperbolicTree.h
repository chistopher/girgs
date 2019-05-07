#pragma once

#include <vector>
#include <utility>
#include <random>
#include <atomic>
#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <mutex>

#include <omp.h>

#include <hypergirgs/ScopedTimer.h>
#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>
#include <hypergirgs/Point.h>
#include <hypergirgs/DistanceFilter.h>
#include <hypergirgs/Generator.h>


namespace hypergirgs {

struct TaskDescription {
    unsigned int cellA;
    unsigned int cellB;
    TaskDescription(unsigned int A, unsigned int B)
        : cellA(A), cellB(B)
    {}
};


template <typename EdgeCallback>
class HyperbolicTree
{
public:

    HyperbolicTree(const std::vector<double>& radii, const std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false);

    void generate(int seed) const;

protected:
    /// Create a set of tasks to be executed in parallel; We'll skip all sampling steps during recursion (call visitCellPairSample!)
    void visitCellPairCreateTasks(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int first_parallel_level,
                                  std::vector<TaskDescription>& parallel_calls) const;

    /// Performs same recursion as visitCellPairCreateTasks, but samples for cells skipp by visitCellPairCreateTasks.
    int visitCellPairSample(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int first_parallel_level,
                                  int num_threads, int thread_shift, default_random_engine& gen) const;

    /// Recursively sample cellA and cellB for level and higher
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, default_random_engine& gen) const;

    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen) const;
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen) const;

    /// takes lower bound on radius for two layers
    unsigned int partitioningBaseLevel(double r1, double r2) const;

    /// 1.0 / connection probability with respect to hyperbolic distance
    double connectionProbRec(double dist) const;

    std::vector<default_random_engine> initialize_prngs(size_t n, unsigned seed) const;

protected:
    EdgeCallback& m_edgeCallback;
    const bool m_profile;

    const size_t m_n; ///< number of nodes

    const double m_coshR; ///< = cosh(R)

    const double m_T; ///< temperature
    const double m_R; ///< radius

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::vector<Point>          m_points;        ///< points ordered by layer first and cell second
    std::vector<unsigned int>   m_first_in_cell; ///< prefix sums into points array
    std::vector<RadiusLayer>    m_radius_layers; ///< data structure to access the points

    std::vector<std::vector<std::pair<unsigned int, unsigned int> > > m_layer_pairs;

    constexpr static size_t filter_size = 100;
    DistanceFilter<filter_size> m_typeI_filter;
    /// filter for layer ij on level l is in  m_typeII_filter[i*m_layers+j][l-2]; -2 because level 0 and 1 have no type 2 cell pairs
    std::vector<std::vector<std::pair<DistanceFilter<filter_size>,DistanceFilter<filter_size>>>> m_typeII_filter;

#ifndef NDEBUG
    mutable long long m_type1_checks{0}; ///< number of node pairs per thread that are checked via a type 1 check
    mutable long long m_type2_checks{0}; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};

template <typename EdgeCallback>
inline HyperbolicTree<EdgeCallback> makeHyperbolicTree(const std::vector<double>& radii, const std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false) {
    return {radii, angles, T, R, edgeCallback, profile};
}

} // namespace hypergirgs

#include <hypergirgs/HyperbolicTree.inl>
