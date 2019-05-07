
#include <hypergirgs/RadiusLayer.h>

#include <cassert>
#include <omp.h>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/ScopedTimer.h>
#include <hypergirgs/IntSort.h>


namespace hypergirgs {

RadiusLayer::RadiusLayer(double r_min, double r_max, unsigned int targetLevel,
                         const Point* base,
                         const unsigned int* prefix_sum)
    : m_r_min{r_min}, 
      m_r_max{r_max}, 
      m_target_level{targetLevel}, 
      m_base{base}, 
      m_prefix_sums{prefix_sum}
{
#ifndef NDEBUG
    const auto cellsInLevel = AngleHelper::numCellsInLevel(targetLevel);

    const auto begin = m_base + m_prefix_sums[0];
    const auto end   = m_base + m_prefix_sums[cellsInLevel];

    for(auto it = begin; it != end; ++it) {
        // check radius lies within layer's radial bound
        assert(m_r_min <= it->radius);
        assert(it->radius <= m_r_max);
    }
#endif
}

std::vector<RadiusLayer> RadiusLayer::buildPartition(const std::vector<double>& radii, const std::vector<double>& angles,
                            const double R, const double layer_height, 
                            std::vector<Point>& points, std::vector<unsigned int>& first_in_cell, // output parameter
                            bool enable_profiling) {

    assert(radii.size() == angles.size());
    assert(layer_height <= R);

    const auto n = radii.size();

    // translate radius <-> layer
    auto radius_to_layer = [R, layer_height] (double radius) {return static_cast<unsigned int>( (R - radius) / layer_height );};
    auto layer_rad_max   = [R, layer_height] (int l) {return R - l * layer_height;};
    auto layer_rad_min   = [R, layer_height] (int l) {return R - l * layer_height - layer_height;};
    const auto r_min_outer = R - layer_height;
    auto num_layers = static_cast<unsigned int>(std::ceil( R / layer_height));

    // generate look-up to get the level of a layer
    const auto level_of_layer = [&] {
        std::vector<int> level_of_layer(num_layers);
        for (int l = 0; l < num_layers; ++l) {
            level_of_layer[l] = partitioningBaseLevel(layer_rad_min(l), r_min_outer, R);
        }
        assert(std::is_sorted(level_of_layer.crbegin(), level_of_layer.crend()));
        return level_of_layer;
    }();

    // since there can be multiple layers at the same level, we cannot
    // rely on AngleHelper::firstCellInLevel to find a unique first cell
    // of a layer. Hence we precompute the first cell as a (reverse)
    // prefix sum
    const auto first_cell_of_layer = [&] {
        std::vector<unsigned int> first_cell_of_layer(num_layers);
        unsigned int sum = 0;
        for (auto l = num_layers; l--;) {
            first_cell_of_layer[l] = sum;
            sum += AngleHelper::numCellsInLevel(level_of_layer[l]);
        }
        return first_cell_of_layer;
    }();
    const auto max_cell_id = first_cell_of_layer.front() + AngleHelper::numCellsInLevel(level_of_layer[0]);

    points = std::vector<Point>(n);
    // pre-compute values for fast distance computation and also compute
    // the cell a point belongs to
    {
        ScopedTimer timer("Classify points & precompute coordinates", enable_profiling);

        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            assert(0 <= radii[i] && radii[i] < R);
            assert(0 <= angles[i] && angles[i] < 2*PI);

            const auto layer = radius_to_layer(radii[i]);
            const auto level = level_of_layer[layer];
            const auto cell = first_cell_of_layer[layer] + AngleHelper::cellForPoint(angles[i], level);
            points[i] = Point(i, radii[i], angles[i], cell);
        }
    }

    // Sort points by cell-ids
    {
        ScopedTimer timer("Sort points", enable_profiling);

        auto compare = [](const Point &a, const Point &b) { return a.cell_id < b.cell_id; };

        intsort::intsort(points, [](const Point &p) { return p.cell_id; }, max_cell_id + 1);
        //alternatively: std::sort(points.begin(), points.end(), compare);

        assert(std::is_sorted(points.begin(), points.end(), compare));
    }

    // prune of empty layers at the back
    for (num_layers = 1; first_cell_of_layer[num_layers - 1] > points[0].cell_id; ++num_layers) {}

    // compute pointers (prefix sums) into points
    constexpr auto gap_cell_indicator = std::numeric_limits<unsigned int>::max();
    first_in_cell = std::vector<unsigned int>(max_cell_id + 1, gap_cell_indicator);
    {
        ScopedTimer timer("Find first point in cell", enable_profiling);

        // First, we mark the begin of cells that actually contain points
        // and repair the gaps (i.e., empty cells) later. In the mean time,
        // the values of those gaps will remain at gap_cell_indicator.
        first_in_cell[points[0].cell_id] = 0;
        first_in_cell[max_cell_id] = n;

        #pragma omp parallel for
        for (int i = 1; i < n; ++i) {
            if (points[i - 1].cell_id != points[i].cell_id) {
                first_in_cell[points[i].cell_id] = i;
            }
        }

        // Now repair gaps: since first_in_cell shall contain
        // a prefix sum, we simply replace any "gap_cell_indicator"
        // with its nearest non-gap successor. In the main loop,
        // this is always the direct successors since we're iterating
        // from right to left.
        #pragma omp parallel
        {
            const auto threads = omp_get_num_threads();
            const auto rank = omp_get_thread_num();
            const auto chunk_size = (max_cell_id + threads - 1) / threads; // = ceil(max_cell_id / threads)

            // Fix right-most element (if gap) of each thread's elements by looking into chunk of next thread.
            // We do not need an end of array check, since it's guaranteed that the last element is n.
            // We're using on this very short code block to avoid UB even if we're only performing word-wise updates.
            #pragma omp single
            {
                for (int r = 0; r < threads; r++) {
                    const auto end = std::min(max_cell_id, chunk_size * (r + 1));
                    int first_non_invalid = end - 1;
                    while (first_in_cell[first_non_invalid] == gap_cell_indicator)
                        first_non_invalid++;
                    first_in_cell[end - 1] = first_in_cell[first_non_invalid];
                }
            }

            const auto begin = std::min(max_cell_id, chunk_size * rank);

            auto i = std::min(max_cell_id, begin + chunk_size);
            while (i-- > begin) {
                first_in_cell[i] = std::min(
                    first_in_cell[i],
                    first_in_cell[i + 1]);
            }
        }

        #ifndef NDEBUG
        {
            assert(points[n-1].cell_id < max_cell_id);

            // assert that we have a prefix sum starting at 0 and ending in n
            assert(first_in_cell[0] == 0);
            assert(first_in_cell[max_cell_id] == n);
            assert(std::is_sorted(first_in_cell.begin(), first_in_cell.end()));

            // check that each point is in its right cell (and that the cell boundaries are correct)
            for (auto cid = 0u; cid != max_cell_id; ++cid) {
                const auto begin = first_in_cell[cid];
                const auto end   = first_in_cell[cid + 1];
                for (auto idx = begin; idx != end; ++idx)
                    assert(points[idx].cell_id == cid);
            }
        }
        #endif
    }

    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    std::vector<RadiusLayer> radius_layers;
    radius_layers.reserve(num_layers);
    {
        ScopedTimer timer("Build data structure", enable_profiling);
        for (auto layer = 0u; layer < num_layers; ++layer) {
            radius_layers.emplace_back(
                layer_rad_min(layer), layer_rad_max(layer), level_of_layer[layer], 
                points.data(), first_in_cell.data() + first_cell_of_layer[layer]);
        }
    }

    return radius_layers;
}

} // namespace hypergirgs
