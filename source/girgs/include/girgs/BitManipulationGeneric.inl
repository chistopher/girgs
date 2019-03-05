#pragma once
#include <array>
#include <cassert>
#include <cstdint>

#ifdef USE_BMI2
#include <immintrin.h>
#endif    


namespace girgs {
namespace BitManipulationDetails {
namespace Generic {

// Generic Implementation of Extractting
template <size_t D>
struct Extract {
    static std::array<uint32_t, D> extract(uint32_t cell) {
        std::array<uint32_t, D> result;

        if (D == 1) {
            result.front() = cell;
            return result;
        }

        std::fill_n(result.begin(), D, 0);

        for(int j = 0; j < (31 + D) / D; j++) {
            for(auto d = 0u; d != D; d++, cell >>= 1) {
                result[d] |= (cell & 1) << j;
            }
        }

        return result;
    }
};

template <>
struct Extract<2> {
    static std::array<uint32_t, 2> extract(uint32_t cell) {
        auto x = cell | (static_cast<uint64_t>(cell & ~1) << 31);

        x = x & 0x5555555555555555;
        x = (x | (x >> 1)) & 0x3333333333333333;
        x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0f;
        x = (x | (x >> 4)) & 0x00ff00ff00ff00ff;
        x = (x | (x >> 8)) & 0x0000ffff0000ffff;

        std::array<uint32_t, 2> result;
        result[0] = static_cast<uint32_t>(x);
        result[1] = static_cast<uint32_t>(x >> 32);

        return result;
    }
};

template <>
struct Extract<4> {
    static std::array<uint32_t, 4> extract(uint32_t cell) {
        auto x = (cell >> 0) | (static_cast<uint64_t>(cell & ~1) << 31);
        auto y = (cell >> 2) | (static_cast<uint64_t>(cell & ~7) << 29);

        x = x & 0x1111111111111111;
        y = y & 0x1111111111111111;
        x = (x | (x >>  3)) & 0x0303030303030303;
        y = (y | (y >>  3)) & 0x0303030303030303;
        x = (x | (x >>  6)) & 0x000f000f000f000f;
        y = (y | (y >>  6)) & 0x000f000f000f000f;
        x = (x | (x >> 12)) & 0x000000ff000000ff;
        y = (y | (y >> 12)) & 0x000000ff000000ff;

        std::array<uint32_t, 4> result;
        result[0] = static_cast<uint32_t>(x);
        result[1] = static_cast<uint32_t>(x >> 32);
        result[2] = static_cast<uint32_t>(y);
        result[3] = static_cast<uint32_t>(y >> 32);
        return result;
    }
};

// Generic Implementation of Interleaving (special case for certain D below)
template <size_t D>
struct Deposit {
    static uint32_t deposit(const std::array<uint32_t, D>& coords) {
        if (D == 1)
            return coords.front();

        unsigned int result = 0u;
        unsigned int bit = 0;

        for(auto l = 0u; l*D < 32 + D; l++) {
            for(auto d = 0u; d != D; d++) {
                result |= ((coords[d] >> l) & 1) << bit++;
            }
        }

        return result;
    }
};

template <>
struct Deposit<2> {
    static uint32_t deposit(const std::array<uint32_t, 2>& coords) {
#ifndef NDEBUG
        for(auto x : coords) assert(x <= 0xffff);
#endif

        uint64_t z = coords[0] | (static_cast<uint64_t>(coords[1]) << 32);

        z = (z | (z << 8)) & 0x00FF00FF00FF00FF;
        z = (z | (z << 4)) & 0x0F0F0F0F0F0F0F0F;
        z = (z | (z << 2)) & 0x3333333333333333;
        z = (z | (z << 1)) & 0x5555555555555555;

        z = z | (z >> 31);

        return static_cast<uint32_t>(z);
    }
};


template <>
struct Deposit<4> {
    static uint32_t deposit(const std::array<uint32_t, 4>& coords) {
#ifndef NDEBUG
        for(auto x : coords) assert(x <= 0xff);
#endif

        uint64_t x = coords[0] | (static_cast<uint64_t>(coords[1]) << 32);
        uint64_t y = coords[2] | (static_cast<uint64_t>(coords[3]) << 32);

        x = (x | (x << 12)) & 0x000F000F000F000F;
        y = (y | (y << 12)) & 0x000F000F000F000F;
        x = (x | (x <<  6)) & 0x0303030303030303;
        y = (y | (y <<  6)) & 0x0303030303030303;
        x = (x | (x <<  3)) & 0x1111111111111111;
        y = (y | (y <<  3)) & 0x1111111111111111;

        x |= x >> 31;
        y |= y >> 31;

        return static_cast<uint32_t>(x | (y << 2));
    }
};

template <unsigned D>
struct Implementation {
    static constexpr unsigned kDimensions = D;

    static uint32_t deposit(const std::array<uint32_t, D>& coords) {
        return Deposit<D>::deposit(coords);
    }

    static std::array<uint32_t, kDimensions> extract(uint32_t cell) {
        return Extract<D>::extract(cell);
    }

    static std::string name() {
        return "Generic";
    }
};

} // namespace Generic
} // namespace BitManipulationDetails
} // namespace girgs
