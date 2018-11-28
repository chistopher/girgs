
#pragma once

#include <vector>
#include <cmath>

#include <girgs/Node.h>
#include <girgs/SpatialTreeCoordinateHelper.h>


namespace girgs {


/**
 * @brief
 *  This class implements the data structure to manage point access described in the paper (Lemma 4.1).
 *  The partitioning of the ground space (Lemma 4.2) implicitly results from the implementation of SpatialTree.
 *
 * @tparam D
 *  the dimension of the geometry
 */
template<unsigned int D>
class WeightLayer {
public:

    WeightLayer() = delete;
    WeightLayer(unsigned int layer, unsigned int targetLevel, const SpatialTreeCoordinateHelper<D>& helper, std::vector<Node*>&& nodes);


    /**
     * @brief
     *  Returns the number of points of this weight layer in a cell.
     *  In the notation of the paper, this function returns \f$ |V_i^{cell}| \f$, where i is the index of this weight layer.
     *
     * @param cell
     *  The cell that contains the points.
     * @param level
     *  The level of the given cell. This should be less or equal to the target level of this weight layer.
     * @return
     *  Returns how many points there are in cells {begin..end} using prefix sums. Begin and end are the first/last descendants of cell in target level.
     */
    int pointsInCell(unsigned int cell, unsigned int level) const;


    /**
     * @brief
     *  Implements the second operation required for the data structure in the paper (Lemma 4.1).
     *  The method finds the first descendant of the given cell in the target level.
     *  Then #m_prefix_sums and #m_A are used to find the requested point.
     *
     * @param cell
     *  The cell that contains the points.
     * @param level
     *  The level of the given cell.
     *  This should be less or equal to the target level of this weight layer.
     * @param k
     *  The point we want to access.
     *  This should be less than the number of points of this weight layer in the given cell
     *  (i.e. less than what was returned by pointsInCell(unsigned int, unsigned int) const ).
     * @return
     *  Returns a pointer to the requested node.
     */
    Node* kthPoint(unsigned int cell, unsigned int level, int k) const;

protected:

    const unsigned int m_layer;             ///< the index of the layer
    const unsigned int m_target_level;      ///< the insertion level for the current weight layer (v(i) = wiw0/W)

    std::vector<Node*> m_nodes;             ///< all nodes of the current weight layer

    std::vector<int>   m_points_in_cell;    ///< the number of points in each cell of target_level
    std::vector<int>   m_prefix_sums;       ///< for each cell c in target level: the sum of points of this layer in all cells <c
    std::vector<Node*> m_A;                 ///< m_A[m_prefix_sums[i]+k] contains the k-th point in the i-th cell of target level
};



} // namespace girgs

#include <girgs/WeightLayer.inl>
