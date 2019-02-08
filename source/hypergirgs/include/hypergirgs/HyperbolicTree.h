#pragma once

#include <vector>
#include <utility>
#include <random>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>
#include <hypergirgs/Point.h>

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

    HyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false);

    void generate(int seed);

protected:
    /// Create a set of tasks to be executed in parallel; We'll skip all sampling steps during recursion (call visitCellPairSample!)
    void visitCellPairCreateTasks(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int first_parallel_level,
                                  std::vector<TaskDescription>& parallel_calls);

    /// Performs same recursion as visitCellPairCreateTasks, but samples for cells skipp by visitCellPairCreateTasks.
    int visitCellPairSample(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int first_parallel_level,
                                  int num_threads, int thread_shift, default_random_engine& gen);

    /// Recursively sample cellA and cellB for level and higher
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, default_random_engine& gen);

    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen);
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen);

    /// takes lower bound on radius for two layers
    unsigned int partitioningBaseLevel(double r1, double r2);

    /// 1.0 / connection probability with respect to hyperbolic distance
    double connectionProbRec(double dist) const;

    /// invConnectionProb(1.0/connectionProbRec(x)) = x
    double invConnectionProb(double p) const;

    template<size_t kFilterStages>
    std::pair<std::array<double, kFilterStages+1>, double> computeFilterStages(double maxProb) const;

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

    std::vector<default_random_engine> initialize_prngs(size_t n, unsigned seed) const;


#ifndef NDEBUG
    long long m_type1_checks{0}; ///< number of node pairs per thread that are checked via a type 1 check
    long long m_type2_checks{0}; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};

template <typename EdgeCallback>
inline HyperbolicTree<EdgeCallback> makeHyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false) {
    return {radii, angles, T, R, edgeCallback, profile};
}

} // namespace hypergirgs

#include <hypergirgs/HyperbolicTree.inl>
