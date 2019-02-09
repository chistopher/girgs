#include <atomic>
#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <mutex>
#include <omp.h>


#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/ScopedTimer.h>

namespace hypergirgs {

template <typename EdgeCallback>
HyperbolicTree<EdgeCallback>::HyperbolicTree(std::vector<double> &radii, std::vector<double> &angles,
    double T, double R, EdgeCallback& edgeCallback, bool enable_profiling)
    : m_edgeCallback(edgeCallback)
    , m_profile(enable_profiling)
    , m_n(radii.size())
    , m_coshR(std::cosh(R))
    , m_T(T)
    , m_R(R)
{
    const auto layer_height = 1.0;

    // compute partition and transfer into own object
    m_radius_layers = RadiusLayer::buildPartition(radii, angles, R, layer_height, enable_profiling);
    m_layers = m_radius_layers.size();
    m_levels = m_radius_layers[0].m_target_level + 1;

    // determine which layer pairs to sample in which level
    {
        ScopedTimer timer("Layer Pairs", enable_profiling);
        m_layer_pairs.resize(m_levels);
        for (auto i = 0u; i < m_layers; ++i)
            for (auto j = 0u; j < m_layers; ++j)
                m_layer_pairs[partitioningBaseLevel(m_radius_layers[i].m_r_min, m_radius_layers[j].m_r_min)].emplace_back(i, j);
    }
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::generate(int seed) {
    const auto num_threads = omp_get_max_threads();
    if(num_threads == 1) {
        default_random_engine master_gen(seed >= 0 ? seed : std::random_device{}());
        visitCellPair(0,0,0, master_gen);
        assert(m_type1_checks + m_type2_checks == static_cast<long long>(m_n-1) * m_n);
        return;
    }

    // We select the first parallel level, s.t. each thread gets on average two tasks
    // with matching cellA==cellB. Those tasks are much more expensive than tasks with cellA != cellB.
    const auto first_parallel_level = static_cast<int>(ceil(log2(2*num_threads)));
    assert(first_parallel_level >= 2);
    const auto num_tasks = (first_parallel_level == 2)
        ? 10
        : ((2 << first_parallel_level) + (3 << (first_parallel_level - 1)));
    std::cout << "First Parallel Level: " << first_parallel_level << "\n";

    // prepare seed_seq and initialize a gen per thread for initial sampling
    auto gens = initialize_prngs(num_threads - 1 + num_tasks,
        seed >= 0 ? seed : std::random_device{}());

    // We have to implement our own task queue, here's the state:
    std::vector<TaskDescription> tasks;

    // allow processing of pq
    std::atomic<bool> tasks_generated{false};
    std::mutex mutex;
    std::condition_variable cv;

    std::atomic<unsigned int> next_task_to_process{0};

    #pragma omp parallel num_threads(num_threads)
    {

        const auto tid = omp_get_thread_num();

    // 1. Phase
        // one thread will generate the task list
        if (tid + 1 == num_threads && first_parallel_level < m_levels) {
            ScopedTimer timer("Gen Tasks");
            tasks.reserve(num_tasks);

            visitCellPairCreateTasks(0, 0, 0, first_parallel_level, tasks);

            // In spirit of the LPT scheduling, we place expensive tasks to be processed first
            std::partition(tasks.begin(), tasks.end(), [] (const TaskDescription& t) {
                return t.cellA == t.cellB;});

            std::cout << "Num Parallel Tasks: " << tasks.size() << "\n";
            assert(num_tasks == tasks.size());

            // inform consumers that tasks are ready
            tasks_generated = true;
            cv.notify_all();
        }

        // all others will sample the cells in the first levels of the recursion tree
        if (tid + 1 < num_threads) {
            visitCellPairSample(0, 0, 0, first_parallel_level, num_threads - 1, tid, gens[tid]);

            // wait until tasks are ready
            if (!tasks_generated) {
                std::unique_lock<std::mutex> lock(mutex);
                cv.wait(lock, [&] {return tasks_generated.load();});
            }
        }

    // 2. Phase
        // we're not using omp for to ensure that the first tasks are processed first
        while(true) {
            const auto i = next_task_to_process.fetch_add(1);
            if (i >= tasks.size()) break;

            auto &task = tasks[i];
            visitCellPair(task.cellA, task.cellB, first_parallel_level, gens[num_threads - 1 + i]);
        }
    }

    assert(m_type1_checks + m_type2_checks == static_cast<long long>(m_n-1) * m_n);
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, default_random_engine& gen) {

    if(!AngleHelper::touching(cellA, cellB, level))
    {   // not touching cells
        // sample all type 2 occurrences with this cell pair
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second, gen);
        return;
    }

    // touching cells

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second, gen);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    auto fA = AngleHelper::firstChild(cellA);
    auto fB = AngleHelper::firstChild(cellB);
    visitCellPair(fA + 0, fB + 0, level+1, gen);
    visitCellPair(fA + 0, fB + 1, level+1, gen);
    visitCellPair(fA + 1, fB + 1, level+1, gen);
    if(cellA != cellB)
        visitCellPair(fA + 1, fB + 0, level+1, gen); // if A==B we already did this call 3 lines above
}

template<typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::visitCellPairCreateTasks(unsigned int cellA, unsigned int cellB,
                                                             unsigned int level,
                                                             unsigned int first_parallel_level,
                                                             std::vector<TaskDescription>& parallel_calls) {

    if(!AngleHelper::touching(cellA, cellB, level))
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    auto fA = AngleHelper::firstChild(cellA);
    auto fB = AngleHelper::firstChild(cellB);

    if(level+1 != first_parallel_level) {
        visitCellPairCreateTasks(fA + 0, fB + 0, level + 1, first_parallel_level, parallel_calls);
        visitCellPairCreateTasks(fA + 0, fB + 1, level + 1, first_parallel_level, parallel_calls);
        visitCellPairCreateTasks(fA + 1, fB + 1, level + 1, first_parallel_level, parallel_calls);
        if (cellA != cellB)
            visitCellPairCreateTasks(fA + 1, fB + 0, level + 1, first_parallel_level, parallel_calls); // if A==B we already did this call 3 lines above
    } else {
        // store tasks
        parallel_calls.emplace_back(fA+0, fB+0);
        parallel_calls.emplace_back(fA+0, fB+1);
        parallel_calls.emplace_back(fA+1, fB+1);
        if (cellA != cellB)
            parallel_calls.emplace_back(fA+1, fB);
    }
}

template<typename EdgeCallback>
int HyperbolicTree<EdgeCallback>::visitCellPairSample(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int first_parallel_level,
                                                                int num_threads, int thread_shift, default_random_engine& gen) {

    auto isMyTurn = [&] {
        if (++thread_shift == num_threads) {
            thread_shift = 0;
            return true;
        }
        return false;
    };

    if(!AngleHelper::touching(cellA, cellB, level))
    {   // not touching cells
        // sample all type 2 occurrences with this cell pair
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                if (isMyTurn())
                    sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second, gen);

        return thread_shift;
    }

    // touching cells

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            if (isMyTurn())
                sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second, gen);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return thread_shift;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    if(level+1 != first_parallel_level) {
        auto fA = AngleHelper::firstChild(cellA);
        auto fB = AngleHelper::firstChild(cellB);
        thread_shift = visitCellPairSample(fA + 0, fB + 0, level + 1, first_parallel_level, num_threads, thread_shift, gen);
        thread_shift = visitCellPairSample(fA + 0, fB + 1, level + 1, first_parallel_level, num_threads, thread_shift, gen);
        thread_shift = visitCellPairSample(fA + 1, fB + 1, level + 1, first_parallel_level, num_threads, thread_shift, gen);
        if (cellA != cellB)
            thread_shift = visitCellPairSample(fA + 1, fB + 0, level + 1, first_parallel_level, num_threads, thread_shift, gen);
    }

    return thread_shift;
}


template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen) {
    auto rangeA = m_radius_layers[i].cellIterators(cellA, level);
    auto rangeB = m_radius_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;

#ifndef NDEBUG
    {
        const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
        const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);
        #pragma omp atomic
        m_type1_checks += (cellA == cellB && i == j) ? sizeV_i_A * (sizeV_i_A - 1)  // all pairs in AxA without {v,v}
                                                     : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
    }
#endif // NDEBUG

    const auto threadId = omp_get_thread_num();

    int kA = 0;
    std::uniform_real_distribution<> dist;
    const auto filters = computeFilterStages<3>(1.0);

    for(auto pointerA = rangeA.first; pointerA != rangeA.second; ++kA, ++pointerA) {
        auto offset = (cellA == cellB && i==j) ? kA+1 : 0;
        for (auto pointerB = rangeB.first + offset; pointerB != rangeB.second; ++pointerB) {
            const auto& nodeInA = *pointerA;
            const auto& nodeInB = *pointerB;

            // pointer magic gives same results
            assert(nodeInA == m_radius_layers[i].kthPoint(cellA, level, kA));
            assert(nodeInB == m_radius_layers[j].kthPoint(cellB, level, std::distance(rangeB.first, pointerB) ));

            // points are in correct cells
            assert(cellA - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInA.angle, level));
            assert(cellB - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInB.angle, level));

            // points are in correct radius layer
            assert(m_radius_layers[i].m_r_min < nodeInA.radius && nodeInA.radius <= m_radius_layers[i].m_r_max);
            assert(m_radius_layers[j].m_r_min < nodeInB.radius && nodeInB.radius <= m_radius_layers[j].m_r_max);

            assert(nodeInA != nodeInB);
            if(m_T==0) {
                if (nodeInA.isDistanceBelowR(nodeInB, m_coshR)) {
                    assert(hyperbolicDistance(nodeInA.radius, nodeInA.angle, nodeInB.radius, nodeInB.angle) < m_R);
                    m_edgeCallback(nodeInA.id, nodeInB.id, threadId);
                }
            } else {
                const auto rnd = dist(gen);

                auto real_dist_cosh = nodeInA.hyperbolicDistanceCosh(nodeInB);
                if (real_dist_cosh > filters.first[static_cast<int>(filters.second * rnd)]) {
                    assert(rnd * connectionProbRec(std::acosh(real_dist_cosh)) >= 1.0);
                    continue;
                }

                if(rnd * connectionProbRec(std::acosh(real_dist_cosh)) < 1.0) {
                    m_edgeCallback(nodeInA.id, nodeInB.id, threadId);
                }
            }
        }
    }
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j, default_random_engine& gen) {

    // TODO use cell iterators
    const auto sizeV_i_A = static_cast<long long>(m_radius_layers[i].pointsInCell(cellA, level));
    const auto sizeV_j_B = static_cast<long long>(m_radius_layers[j].pointsInCell(cellB, level));
    if (m_T == 0 || sizeV_i_A == 0 || sizeV_j_B == 0) {
#ifndef NDEBUG
        #pragma omp atomic
        m_type2_checks += 2ll * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG
        return;
    }

    // get upper bound for probability
    auto r_boundA = m_radius_layers[i].m_r_min;
    auto r_boundB = m_radius_layers[j].m_r_min;
    auto angular_distance_lower_bound = AngleHelper::dist(cellA, cellB, level);
    auto dist_lower_bound = hyperbolicDistance(r_boundA, 0, r_boundB, angular_distance_lower_bound);
    auto max_connection_prob = 1.0 / connectionProbRec(dist_lower_bound);

    // if we must sample all pairs we treat this as type 1 sampling
    // also, 1.0 is no valid prob for a geometric dist (see c++ std)
    if(max_connection_prob == 1.0){
        sampleTypeI(cellA, cellB, level, i, j, gen);
        return;
    }

    const auto num_pairs = sizeV_i_A * sizeV_j_B;
    const auto expected_samples = num_pairs * max_connection_prob;

#ifndef NDEBUG
    #pragma omp atomic
    m_type2_checks += 2llu * num_pairs;
#endif // NDEBUG

    if(expected_samples < 1e-6)
        return;

    // init geometric distribution
    const auto threadId = omp_get_thread_num();
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);
    std::uniform_real_distribution<> dist(0.0, max_connection_prob);

    if (expected_samples > 10.) {
        const auto filters = computeFilterStages<3>(max_connection_prob);

        for (auto r = geo(gen); r < num_pairs; r += 1 + geo(gen)) {
            // determine the r-th pair
            const auto& nodeInA = m_radius_layers[i].kthPoint(cellA, level, r%sizeV_i_A);
            const auto& nodeInB = m_radius_layers[j].kthPoint(cellB, level, r/sizeV_i_A);

            // points are in correct cells
            assert(cellA - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInA.angle, level));
            assert(cellB - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInB.angle, level));

            // points are in correct radius layer
            assert(m_radius_layers[i].m_r_min < nodeInA.radius && nodeInA.radius <= m_radius_layers[i].m_r_max);
            assert(m_radius_layers[j].m_r_min < nodeInB.radius && nodeInB.radius <= m_radius_layers[j].m_r_max);

            // get actual connection probability
            const auto real_dist_cosh = nodeInA.hyperbolicDistanceCosh(nodeInB);
            assert(angular_distance_lower_bound <= std::abs(nodeInA.angle - nodeInB.angle));
            assert(angular_distance_lower_bound <= std::abs(nodeInB.angle - nodeInA.angle));
            assert(std::acosh(real_dist_cosh) >= dist_lower_bound);
            assert(std::acosh(real_dist_cosh) > m_R);

            const auto rnd = dist(gen);
            if (real_dist_cosh > filters.first[static_cast<int>(filters.second * rnd)]) {
                assert(rnd * connectionProbRec(std::acosh(real_dist_cosh)) >= 1.0);
                continue;
            }

            auto connection_prob = connectionProbRec(std::acosh(real_dist_cosh));
            if(rnd * connection_prob < 1.0) {
                m_edgeCallback(nodeInA.id, nodeInB.id, threadId);
            }
        }

    } else {
        for (auto r = geo(gen); r < num_pairs; r += 1 + geo(gen)) {
            // determine the r-th pair
            const auto& nodeInA = m_radius_layers[i].kthPoint(cellA, level, r%sizeV_i_A);
            const auto& nodeInB = m_radius_layers[j].kthPoint(cellB, level, r/sizeV_i_A);

            // points are in correct cells
            assert(cellA - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInA.angle, level));
            assert(cellB - AngleHelper::firstCellOfLevel(level) == AngleHelper::cellForPoint(nodeInB.angle, level));

            // points are in correct radius layer
            assert(m_radius_layers[i].m_r_min < nodeInA.radius && nodeInA.radius <= m_radius_layers[i].m_r_max);
            assert(m_radius_layers[j].m_r_min < nodeInB.radius && nodeInB.radius <= m_radius_layers[j].m_r_max);

            // get actual connection probability
            const auto real_dist = nodeInA.hyperbolicDistance(nodeInB);
            assert(angular_distance_lower_bound <= std::abs(nodeInA.angle - nodeInB.angle));
            assert(angular_distance_lower_bound <= std::abs(nodeInB.angle - nodeInA.angle));
            assert(real_dist >= dist_lower_bound);
            assert(real_dist > m_R);

            const auto connection_prob = connectionProbRec(real_dist);
            if(dist(gen) * connection_prob < 1.0) {
                m_edgeCallback(nodeInA.id, nodeInB.id, threadId);
            }
        }
    }
}


template <typename EdgeCallback>
std::vector<default_random_engine> HyperbolicTree<EdgeCallback>::initialize_prngs(size_t n, unsigned seed) const {
    std::vector<default_random_engine> gens;

    // we need a generator for each tasks and also for each but one threads
    constexpr auto state_size = default_random_engine::state_size;

    // get high quality seed values from one PRNG
    std::vector<uint32_t> seeds(n * state_size);
    {
        default_random_engine master_gen(seed);
        std::uniform_int_distribution<std::uint32_t> dist;
        std::generate(seeds.begin(), seeds.end(), [&] {return dist(master_gen);});
    }

    // use a seed_seq to initialize all prngs
    gens.reserve(n);
    for(int i=0; i < n; ++i) {
        auto begin = std::next(seeds.begin(), i * state_size);
        auto end = std::next(begin, state_size);
        std::seed_seq ss(begin, end);
        gens.emplace_back(ss);
    }

    return gens;
}

template <typename EdgeCallback>
unsigned int HyperbolicTree<EdgeCallback>::partitioningBaseLevel(double r1, double r2) {
    return RadiusLayer::partitioningBaseLevel(r1, r2, m_R);
}

template<typename EdgeCallback>
double HyperbolicTree<EdgeCallback>::connectionProbRec(double dist) const {
    return 1.0 + std::exp(0.5/m_T*(dist-m_R));
}

template<typename EdgeCallback>
double HyperbolicTree<EdgeCallback>::invConnectionProb(double p) const {
    return m_R + 2*m_T*std::log(1.0 / p - 1);
}

template<typename EdgeCallback>
template<size_t kFilterStages>
std::pair<std::array<double, kFilterStages+1>, double> HyperbolicTree<EdgeCallback>::computeFilterStages(double max_connection_prob) const {
    std::array<double, kFilterStages+1> filters;
    for(int i=0; i <= kFilterStages; i++)
        filters[i] = cosh(invConnectionProb(max_connection_prob / kFilterStages * i));
    const double filter_width = kFilterStages / max_connection_prob;
    return {filters, filter_width};
}


} // namespace hypergirgs
