#pragma once

#include <cassert>
#include <algorithm>
#include <cmath>

#ifndef NDEBUG
#define POINT_WITH_ORIGINAL
#endif

namespace hypergirgs {

/*** TODO docs
 *
 * @param r1
 * @param phi1
 * @param r2
 * @param phi2
 * @return
 */
static double hyperbolicDistance(double r1, double phi1, double r2, double phi2) {
    return acosh(std::max(1., cosh(r1 - r2) + (1. - cos(phi1 - phi2)) * sinh(r1) * sinh(r2)));
}

struct Point {
    Point() {}; // prevent initialization of members
    Point(const int id, const double radius, const double angle, int cell_id = 0) :
          id{id}
        , cell_id{cell_id}
        , invsinh_r{1.0 / std::sinh(radius)}
        , coth_r{std::cosh(radius) / std::sinh(radius)}
        , cos_phi{std::cos(angle)}
        , sin_phi{std::sin(angle)}
#ifdef POINT_WITH_ORIGINAL
        , radius{radius}
        , angle{angle}
#endif // POINT_WITH_ORIGINAL
    {
        assert(0 <= angle && angle < 2*3.14159265358979323846);
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

    /// Returns cosh(hyperbolicDistance to pt)
    double hyperbolicDistanceCosh(const Point& pt) const noexcept {
        // cosh(dist)
        // = cosh(r1)cosh(r2) - sinh(r1)sinh(r2)*cos(p1-p2)
        // = cosh(r1)cosh(r2) - sinh(r1)sinh(r2)*(sin(p1)sin(p2)+cos(p1)cos(p2))
        // = cosh(r1)/sinh(r1)cosh(r2)/sinh(r2) * sinh(r1)sinh(r2) - sinh(r1)sinh(r2)*(sin(p1)sin(p2)+cos(p1)cos(p2))
        // = sinh(r1)sinh(r2) * (cosh(r1)/sinh(r1)cosh(r2)/sinh(r2) - (sin(p1)sin(p2)+cos(p1)cos(p2)))
        return std::max(1.0, (coth_r * pt.coth_r - cos_phi * pt.cos_phi - sin_phi * pt.sin_phi) / (invsinh_r * pt.invsinh_r));
    }

    /// Returns hyperbolic distance to pt
    double hyperbolicDistance(const Point& pt) const noexcept {
        return std::acosh(hyperbolicDistanceCosh(pt));
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


#ifdef POINT_WITH_ORIGINAL
    double radius;    ///< = radius
    double angle;     ///< = angle
#endif // POINT_WITH_ORIGINAL

    void prefetch() const noexcept {
#if defined(__GNUC__) || defined(__clang__)
        __builtin_prefetch(&id, 0);
    #ifndef POINT_WITH_ORIGINAL
        __builtin_prefetch(&sin_phi, 0);
    #else
        __builtin_prefetch(&angle, 0);
    #endif
#endif
    }
};

} // namespace hypergirgs
