#include <algorithm>
#include <cassert>
#include <omp.h>

#include <hypergirgs/Hyperbolic.h>
#include <ScopedTimer.h>
#include <IntSort.h>

namespace hypergirgs {

template <typename EdgeCallback>
HyperbolicTree<EdgeCallback>::HyperbolicTree(std::vector<double> &radii, std::vector<double> &angles, double T, double R, EdgeCallback& edgeCallback, bool profile)
    : m_edgeCallback(edgeCallback)
    , m_profile(profile)
    , m_n(radii.size())
    , m_coshR(std::cosh(R))
    , m_T(T)
    , m_R(R)
    , m_points(m_n)
    , m_gen()
    , m_dist()
    #ifndef NDEBUG
    , m_type1_checks(0)
    , m_type2_checks(0)
    #endif // NDEBUG
{
    assert(radii.size() == angles.size());

    // translate radius <-> layer
    const auto layer_height = 1.0; //std::log(2.0)/0.75;
    auto radius_to_layer = [R, layer_height] (double radius) {return static_cast<unsigned int>( (R - radius) / layer_height );};
    auto layer_rad_max   = [R, layer_height] (int l) {return R - l * layer_height;};
    auto layer_rad_min   = [R, layer_height] (int l) {return R - l * layer_height - layer_height;};
    const auto r_min_outer = R - layer_height;

    m_layers = static_cast<unsigned int>(std::ceil(R/layer_height));

    // generate look-up to get the level of a layer
    const auto level_of_layer = [&] {
        std::vector<int> level_of_layer(m_layers);
        for (int l = 0; l != m_layers; ++l) {
            level_of_layer[l] = partitioningBaseLevel(layer_rad_min(l), r_min_outer);
        }
        assert(std::is_sorted(level_of_layer.crbegin(), level_of_layer.crend()));
        return level_of_layer;
    }();

    // since there can be multiple layers at the same level, we cannot
    // rely on AngleHelper::firstCellInLevel find a unique first cell
    // of a layer. Hence we precompute the first cell as a (reverse)
    // prefix sum
    const auto first_cell_of_layer = [&] {
        std::vector<unsigned int> first_cell_of_layer(m_layers);
        unsigned int sum = 0;
        for (auto l = m_layers; l--;) {
            first_cell_of_layer[l] = sum;
            sum += AngleHelper::numCellsInLevel(level_of_layer[l]);
        }
        return first_cell_of_layer;
    }();
    const auto max_cell_id = first_cell_of_layer.front() + AngleHelper::numCellsInLevel(level_of_layer[0]);

    // pre-compute values for fast distance computation and also compute
    // the cell a point belongs to
    {
        ScopedTimer timer("Classify points & precompute coordinates", m_profile);

        #pragma omp parallel for
        for (int i = 0; i < m_n; ++i) {
            const auto layer = radius_to_layer(radii[i]);
            const auto level = level_of_layer[layer];
            const auto cell = first_cell_of_layer[layer] + AngleHelper::cellForPoint(angles[i], level);

            auto &pt = m_points[i];
            pt = Point(i, radii[i], angles[i], cell);
        }
    }

    // Sort points by cell-ids
    {
        ScopedTimer timer("Sort points", m_profile);

        auto compare = [] (const Point& a, const Point& b) {return a.cell_id < b.cell_id;};

        intsort::intsort(m_points, [] (const Point& p) {return p.cell_id;}, max_cell_id + 1);
        //std::sort(m_points.begin(), m_points.end(), compare);

        assert(std::is_sorted(m_points.begin(), m_points.end(), compare));
    }

    // prune of empty layers at the back
    for(m_layers = 1; first_cell_of_layer[m_layers-1] > m_points.front().cell_id; ++m_layers) {}

    // compute pointers into m_points
    {
        ScopedTimer timer("Find first point in cell", m_profile);
        constexpr auto gap_cell_indicator = std::numeric_limits<unsigned int>::max();


        // First we mark the begin of cells that actually contain points
        // and repair the gaps (i.e., empty cells) later.
        // In the mean time, gaps will remain at m_n.
        m_first_point_in_cell[m_points.front().cell_id] = 0;
        for(auto i=1; i != m_n; ++i) {
            if (m_points[i-1].cell_id != m_points[i].cell_id) {
                m_first_point_in_cell[m_points[i].cell_id] = i;
            }
        }

        // Now repair gaps: since m_first_point_in_cell shell contain
        // a prefix sum, we simply replace any m_n (but the last)
        // with its left-most non-m_n successor
        for(auto i=max_cell_id-1; i--;) {
            m_first_point_in_cell[i] = std::min(
                m_first_point_in_cell[i],
                m_first_point_in_cell[i+1]);
        }

#ifndef NDEBUG
        {
            assert(m_points.back().cell_id < max_cell_id);

            // assert that we have a prefix sum starting at 0 and ending in m_n
            assert(m_first_point_in_cell.front() == 0);
            assert(m_first_point_in_cell.back()  == m_n);
            assert(std::is_sorted(m_first_point_in_cell.cbegin(), m_first_point_in_cell.cend()));

            // check that each point is in its right cell (and that the cell boundaries are correct)
            for(auto cid = 0u; cid != max_cell_id; ++cid) {
                const auto begin = m_points.data() + m_first_point_in_cell[cid];
                const auto end   = m_points.data() + m_first_point_in_cell[cid+1];
                for(auto it = begin; it != end; ++it)
                    assert(it->cell_id == cid);
            }
        }
#endif
    }


    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    {
        ScopedTimer timer("Build data structure", m_profile);
        for (auto layer = 0u; layer < m_layers; ++layer) {
            m_radius_layers.emplace_back(
                layer_rad_min(layer), layer_rad_max(layer), level_of_layer[layer],
                m_points.data(), m_first_point_in_cell.data() + first_cell_of_layer[layer]);
        }
    }

    m_levels = m_radius_layers[0].m_target_level + 1;

    // determine which layer pairs to sample in which level
    {
        ScopedTimer timer("Layer Pairs", m_profile);
        m_layer_pairs.resize(m_levels);
        for (auto i = 0u; i < m_layers; ++i)
            for (auto j = 0u; j < m_layers; ++j)
                m_layer_pairs[partitioningBaseLevel(m_radius_layers[i].m_r_min, m_radius_layers[j].m_r_min)].emplace_back(i, j);
    }
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::generate(int seed) {
    m_gen.seed(seed >= 0 ? seed : std::random_device{}());
    m_dist.reset();
    visitCellPair(0,0,0);
    assert(m_type1_checks + m_type2_checks == static_cast<long long>(m_n-1) * m_n);
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {

    if(!AngleHelper::touching(cellA, cellB, level))
    {   // not touching cells
        // sample all type 2 occurrences with this cell pair
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // touching cells

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    auto fA = AngleHelper::firstChild(cellA);
    auto fB = AngleHelper::firstChild(cellB);
    visitCellPair(fA + 0, fB + 0, level+1);
    visitCellPair(fA + 0, fB + 1, level+1);
    visitCellPair(fA + 1, fB + 1, level+1);
    if(cellA != cellB)
        visitCellPair(fA + 1, fB + 0, level+1); // if A==B we already did this call 3 lines above
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j) {
    auto rangeA = m_radius_layers[i].cellIterators(cellA, level);
    auto rangeB = m_radius_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;

#ifndef NDEBUG
    {
        const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
        const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);
        m_type1_checks += (cellA == cellB && i == j) ? sizeV_i_A * (sizeV_i_A - 1)  // all pairs in AxA without {v,v}
                                                     : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
    }
#endif // NDEBUG

    const auto threadId = omp_get_thread_num();

    int kA = 0;
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

            // points are in correct weight layer
            assert(m_radius_layers[i].m_r_min < nodeInA.radius && nodeInA.radius <= m_radius_layers[i].m_r_max);
            assert(m_radius_layers[j].m_r_min < nodeInB.radius && nodeInB.radius <= m_radius_layers[j].m_r_max);

            assert(nodeInA != nodeInB);
            if (nodeInA.isDistanceBelowR(nodeInB, m_coshR)) {
                assert(hyperbolicDistance(nodeInA.radius, nodeInA.angle, nodeInB.radius, nodeInB.angle) < m_R);
                m_edgeCallback(nodeInA.id, nodeInB.id, threadId);
            }
        }
    }
}

template <typename EdgeCallback>
void HyperbolicTree<EdgeCallback>::sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j) {

    auto sizeV_i_A = m_radius_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_radius_layers[j].pointsInCell(cellB, level);
    if (sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

#ifndef NDEBUG
    m_type2_checks += 2llu * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    if (m_T == 0)
        return;

    // TODO add content

}

template <typename EdgeCallback>
unsigned int HyperbolicTree<EdgeCallback>::partitioningBaseLevel(double r1, double r2) {
    auto level = 0u;
    auto cellDiameter = 2.0*PI;
    // find deepest level in which points in all non-touching cells are not connected
    while(hypergirgs::hyperbolicDistance(r1, 0, r2, (cellDiameter/2)) > m_R){
        level++;
        cellDiameter /= 2;
    }
    return level;
}


} // namespace hypergirgs
