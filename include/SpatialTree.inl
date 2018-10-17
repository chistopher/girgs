

#include <memory>


template<unsigned int D, unsigned int L>
constexpr SpatialTree<D,L>::SpatialTree()
    : m_coords()
    , m_coords2Index()
    , m_points_in_cell()
    , m_prefix_sums()
    , m_A(nullptr)
    , m_maxLevel(0)
    , m_weight_layers(nullptr)
{

    // calculate coords // TODO find implicit way to get from index to coords and back
    for(auto d=0; d<D; ++d) m_coords[0][d] = 0; // first cell has [0]^D // no "m_coords[0].fill(0);" because that is not constexpr

    for(auto l=1; l<L; ++l) {
        for(auto cell = firstCellOfLevel(l); cell < firstCellOfLevel(l+1); ++cell) {
            m_coords[cell] = m_coords[parent(cell)];
            auto bitmask = ((cell-1)&(numChildren-1)); // & is unnecessary but removes all but the last D bits
            for(auto d=0; d<D; ++d) { // check orientation to parent in all dimensions with d-th bit in bitmask
                m_coords[cell][d] *= 2; // every level doubles each dimension
                m_coords[cell][d] += static_cast<bool>(bitmask&(1<<d)); // add 1 if we are the "right" of parent in dimension d
            }

            // save inverse mapping
            auto packedCoord = 0;
            for(auto d=0; d<D; ++d)
                packedCoord += m_coords[cell][d] * (1<<(l*d));
            m_coords2Index[firstCellOfLevel(l) + packedCoord] = cell;
        }
    }
}


template<unsigned int D, unsigned int L>
std::array<std::pair<double, double>, D> SpatialTree<D, L>::bounds(int cell, unsigned int level) const {
    auto diameter = 1.0 / (1<<level);
    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0; d<D; ++d)
        result[d]= { m_coords[cell][d]*diameter, (m_coords[cell][d]+1)*diameter };
    return result;
}


template<unsigned int D, unsigned int L>
unsigned int SpatialTree<D, L>::cellForPoint(std::array<double, D>& point, unsigned int targetLevel) const {

    // calculate coords
    auto diameter = 1.0 / (1<<targetLevel);
    std::array<int, D> coords;
    for(auto d=0; d<D; ++d) {
        coords[d] = static_cast<int>(point[d] / diameter);
        assert(coords[d] < (1<<targetLevel));
    }

    // get from coords to index
    auto packedCoord = 0;
    for(auto d=0; d<D; ++d)
        packedCoord += coords[d] * (1<<(targetLevel*d));
    return m_coords2Index[firstCellOfLevel(targetLevel) + packedCoord];
}


template<unsigned int D, unsigned int L>
typename SpatialTree<D, L>::Graph SpatialTree<D, L>::generateGraph(std::vector<int>& weights) {

    auto graph = Graph(weights.size());
    // sample positions
    // TODO

    // build weight layers
    auto weightLayers = std::make_unique(decltype(m_weight_layers)());
    m_weight_layers = weightLayers.get();
    // TODO do stuff
    m_maxLevel = weightLayers.size(); // TODO or -1?

    // insert points in spatial data structure
    auto A = std::make_unique(decltype(m_A)()); // in fact we only need an array of size m_maxLevel <= L
    m_A = A.get();
    // TODO fill m_points_in_cell, m_prefix_sums, m_A
    // TODO counting sort for each level / layer

    // sample all edges
    //visitCellPair(0,0, graph);

    return graph;
}


template<unsigned int D, unsigned int L>
void SpatialTree<D, L>::visitCellPair(unsigned int cellA, unsigned int cellB, SpatialTree::Graph &graph) const {
    /*
    let l be level of cells

    if(touching A,B or A=B)
        // sample all type 1 occurrences with this cell pair
        sample edges for weightlayers i,j with (wi*wj/W) somehow determines l as target level // is i+j = l right ??? -> see "we sample all point-pairs" proof
    else // not touching
        // sample all type 2 occurrences with this cell pair
        // the type 2 occurrences can come from any pair of lower levels i,j>=l
        // assert(parent must be touching) to be sure we are type II pair; otherwise we would not have done this recursive call
        for all weightlayer pairs i,j>=l
            sample edges between V_i^A and V_j^B
            // for two points x in V_i^A and y in V_j^B this can happen only one time, because there is just one pair of cells (A', B') in the hierarchy with x in A' and y in B' where:
            // 1. A' and B' are not touching
            // 2: parent(A') touches parent(B')
            // It is easy to see, that no children of these cells will sample those points again, because the recursion ends here

    if(A=B)
        recursive call for all children pairs of cell A=B

    if(touching A,B and A!=B)
        recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they dont
    */
    // TODO implement
}


template<unsigned int D, unsigned int L>
int SpatialTree<D,L>::pointsInCell(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel) const {
    assert(fromLevel <= targetLevel);
    assert(firstCellOfLevel(fromLevel) <= cell && cell < firstCellOfLevel(fromLevel+1)); // cell is from fromLevel

    // we want the begin-th and end-th cell in level targetLevel to be the first and last descendant of cell in this level
    // we could apply the firstChild function to find the first descendant but this is in O(1)
    auto descendants = numCellsInLevel(targetLevel - fromLevel); // 2^(D*L)
    auto localIndexCell = cell - firstCellOfLevel(fromLevel);
    auto localIndexDescendant = localIndexCell * descendants; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendent
    auto begin = localIndexDescendant + firstCellOfLevel(targetLevel);
    auto end = begin + descendants - 1;

    assert(begin < firstCellOfLevel(targetLevel+1));
    assert(end < firstCellOfLevel(targetLevel+1));

    return m_prefix_sums[end] - m_prefix_sums[begin] + m_points_in_cell[end];
}


template<unsigned int D, unsigned int L>
typename SpatialTree<D, L>::Node* SpatialTree<D, L>::kthPoint(unsigned int cell, unsigned int fromLevel, unsigned int targetLevel, int k) const {
    assert(fromLevel <= targetLevel);
    assert(firstCellOfLevel(fromLevel) <= cell && cell < firstCellOfLevel(fromLevel+1)); // cell is from fromLevel

    // same as in "pointsInCell"
    auto descendants = numCellsInLevel(targetLevel - fromLevel);
    auto localIndexCell = cell - firstCellOfLevel(fromLevel);
    auto localIndexDescendant = localIndexCell * descendants;
    auto begin = localIndexDescendant + firstCellOfLevel(targetLevel);

    // TODO think about this!
    return (*m_A)[targetLevel][m_prefix_sums[begin] + k];
}


