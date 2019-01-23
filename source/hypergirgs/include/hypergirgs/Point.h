#pragma once

#include <cassert>
#include <cmath>

#include <hypergirgs/Hyperbolic.h>

namespace hypergirgs {

struct Point {
    Point() {}; // prevent initialization of members
    Point(const int id, const double radius, const double angle, int cell_id = 0) :
          id{id}
        , cell_id{cell_id}
        , invsinh_r{1.0 / std::sinh(radius)}
        , coth_r{std::cosh(radius) / std::sinh(radius)}
        , cos_phi{std::cos(angle)}
        , sin_phi{std::sin(angle)}
        , radius{radius}
        , angle{angle}
    {
        assert(0 <= angle && angle < 2*PI);
        assert(0 <= radius);
        assert(0 <= id);
    }

    /// Check whether distance between this point and point pt is below the threshold R
    /// without using trigonometric functions. (Useful in the threshold model)
    /// @warning Pass cosh(R) rather than R as second parameter!
    bool isDistanceBelowR(const Point& pt, const double coshR) const noexcept {
        assert(coshR > 1.0 / invsinh_r); // should fire eventually if R rather than cosh(R) is passed
        return cos_phi * pt.cos_phi + sin_phi * pt.sin_phi >
            coth_r * pt.coth_r - coshR * invsinh_r * pt.invsinh_r;
    }

    /// Check whether node ids match
    bool operator==(const Point& o) const noexcept {
        return id == o.id;
    }

    /// Check whether node ids are unequal
    bool operator!=(const Point& o) const noexcept {
        return id != o.id;
    }

    int    id;        ///< node id
    int    cell_id;   ///< id of cell node will stored

    double invsinh_r; ///< = 1.0 / sinh(radius)
    double coth_r;    ///< = coth(radius) = cosh(radius) / sinh(radius)
    double cos_phi;   ///< = cos(angle)
    double sin_phi;   ///< = sin(angle)

    // TODO also use improved distance computations for T>0, thus remove these in release again
    double radius;    ///< = radius
    double angle;     ///< = angle
};

} // namespace hypergirgs
