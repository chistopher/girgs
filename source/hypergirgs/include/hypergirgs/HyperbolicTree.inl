#include <cassert>
#include <omp.h>

#include <hypergirgs/Hyperbolic.h>
#include <hypergirgs/HyperbolicTree.h>

namespace hypergirgs {

template <typename EdgeCallback>
HyperbolicTree<EdgeCallback>::HyperbolicTree(std::vector<double> &radii, std::vector<double> &angles, double T, double R, EdgeCallback& edgeCallback)
: m_edgeCallback(edgeCallback)
, m_n(radii.size())
, m_coshR(std::cosh(R))
, m_T(T)
, m_R(R)
, m_gen()
, m_dist()
#ifndef NDEBUG
, m_type1_checks(0)
, m_type2_checks(0)
#endif // NDEBUG
{
    assert(radii.size() == angles.size());

    // pre-compute values for distance
    std::vector<Point> pre_points(radii.size());
    #pragma omp parallel for if (m_n > 10000)
    for (int i = 0; i < radii.size(); ++i) {
        pre_points[i] = Point(i, radii[i], angles[i]);
    }

    // create layer
    m_layers = static_cast<unsigned int>(std::ceil(R));
    auto weightLayerNodes = std::vector<std::vector<int>>(m_layers);
    for(auto i = 0u; i < radii.size(); ++i) // layer i has nodes from (R-i-1 to R-i]
        weightLayerNodes[static_cast<unsigned int>(R-radii[i])].push_back(i);

    // ignore empty layers of higher radius
    for(;m_layers>0;m_layers--)
        if(!weightLayerNodes[m_layers-1].empty())
            break;

    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    for (auto layer = 0u; layer < m_layers; ++layer)
        m_radius_layers.emplace_back(R - layer - 1, R - layer, partitioningBaseLevel(R - layer - 1, R - 1),
                                     std::move(weightLayerNodes[layer]), angles, pre_points);
    m_levels = m_radius_layers[0].m_target_level + 1;

    // determine which layer pairs to sample in which level
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(m_radius_layers[i].m_r_min, m_radius_layers[j].m_r_min)].emplace_back(i,j);
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

    auto sizeV_i_A = m_radius_layers[i].pointsInCell(cellA, level);
    auto sizeV_j_B = m_radius_layers[j].pointsInCell(cellB, level);
    if (sizeV_i_A == 0 || sizeV_j_B == 0)
        return;

#ifndef NDEBUG
    m_type1_checks += (cellA == cellB && i == j)
        ? sizeV_i_A * (sizeV_i_A-1)  // all pairs in AxA without {v,v}
        : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
#endif // NDEBUG

    const auto * firstA = m_radius_layers[i].firstPointPointer(cellA, level);
    const auto * firstB = m_radius_layers[j].firstPointPointer(cellB, level);
    const auto threadId = omp_get_thread_num();

    for(int kA=0; kA<sizeV_i_A; ++kA){
        for (int kB =(cellA == cellB && i==j ? kA+1 : 0); kB<sizeV_j_B; ++kB) {
            const auto& nodeInA = *(firstA+kA);
            const auto& nodeInB = *(firstB+kB);

            // pointer magic gives same results
            assert(nodeInA == m_radius_layers[i].kthPoint(cellA, level, kA));
            assert(nodeInB == m_radius_layers[j].kthPoint(cellB, level, kB));

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
    auto cellDiameter = 2*M_PI;
    // find deepest level in which points in all non-touching cells are not connected
    while(hypergirgs::hyperbolicDistance(r1, 0, r2, (cellDiameter/2)) > m_R){
        level++;
        cellDiameter /= 2;
    }
    return level;
}


} // namespace hypergirgs
