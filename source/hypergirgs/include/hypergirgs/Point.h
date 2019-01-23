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
#ifndef NDEBUG
        , radius{radius}
        , angle{angle}
#endif // NDEBUG
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

    double hyperbolicDistance(const Point& pt) const noexcept {
        // dist =
        // acosh( cosh(r1)cosh(r2) - sinh(r1)sinh(r2)*cos(p1-p2) )
        // acosh( cosh(r1)cosh(r2) - sinh(r1)sinh(r2)*(sin(p1)sin(p2)+cos(p1)cos(p2)) )
        // acosh( cosh(r1)/sinh(r1)cosh(r2)/sinh(r2) * sinh(r1)sinh(r2) - sinh(r1)sinh(r2)*(sin(p1)sin(p2)+cos(p1)cos(p2)) )
        // acosh( sinh(r1)sinh(r2) * (cosh(r1)/sinh(r1)cosh(r2)/sinh(r2) - (sin(p1)sin(p2)+cos(p1)cos(p2))) )
        return std::acosh(std::max(1.0,
                (coth_r * pt.coth_r - cos_phi * pt.cos_phi - sin_phi * pt.sin_phi) / (invsinh_r * pt.invsinh_r)
        ));
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

#ifndef NDEBUG
    double radius;    ///< = radius
    double angle;     ///< = angle
#endif // NDEBUG
};

} // namespace hypergirgs
