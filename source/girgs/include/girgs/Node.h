
#pragma once

#include <array>
#include <vector>
#include <cassert>


namespace girgs {

template<unsigned int D>
struct Node {
    std::array<double, D>   coord;
    double                  weight;
    int                     index;

    Node() = default;
    Node(const std::vector<double>& _coord, double _weight, int _index)
    : weight(_weight), index(_index) {
        assert(_coord.size()==D);
        for(auto d=0u; d<D; ++d)
            coord[d] = _coord[d];
    }

    double distance(const Node& other) const {
        auto result = 0.0;
        for(auto d=0u; d<D; ++d){
            auto dist = std::abs(coord[d] - other.coord[d]);
            dist = std::min(dist, 1.0-dist);
            result = std::max(result, dist);
        }
        return result;
    }
};


} // namespace girgs
