#pragma once

#include <immintrin.h>

namespace girgs {
namespace BitManipulationDetails {
namespace BMI2 {

template <size_t D>
struct Extract {
    static std::array<uint32_t, D> extract(uint32_t cell) {
        std::array<uint32_t, D> result;

        for(int i = 0; i < D; ++i)
            result[i] = _pext_u32(cell, BitPattern<D, uint32_t>::kEveryDthBit << i);

        return result;
    }
};

struct PImpl {
    constexpr static uint32_t highbits(int x) {
        return static_cast<uint32_t>((1llu << x) - 1);
    }

    static uint32_t impl_pdep(uint32_t a) noexcept {
        return a;
    }

    static uint32_t impl_pdep(uint32_t a, uint32_t b) noexcept {
        auto values = (b << 16) | (a & 0xffff);
        auto deposited = _pdep_u64(values, 0x5555555555555555ull);
        return static_cast<uint32_t>(deposited | (deposited >> 31));
    }

    static uint32_t impl_pdep(uint32_t a, uint32_t b, uint32_t c) noexcept {
        const auto rec = impl_pdep(a, b);

        const auto values = ((c & highbits(10)) << 22) | (rec & highbits(22));
        const auto deposited = _pdep_u64(values, 0x24924924db6db6db); // = 0b00'100100100100100100100100100100'11011011011011011011011011011011

        return (deposited & highbits(32)) | (deposited >> 32);
    }

    static uint32_t impl_pdep(uint32_t a, uint32_t b, uint32_t c, uint32_t d) noexcept {
        return impl_pdep(impl_pdep(a, c), impl_pdep(b, d));
    }

    static uint32_t impl_pdep(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e) noexcept {
        const auto rec = impl_pdep(a,b,c,d);

        const auto values = ((e & highbits(6)) << 26) | (rec & highbits(26));
        const auto deposited = _pdep_u64(values, 0x21084210def7bdef); //  0b00'100001000010000100001000010000'11011110111101111011110111101111);

        return (deposited & highbits(32)) | (deposited >> 32);
    }
};

template <unsigned D>
struct Deposit;

#define ImplSpec(X, ...) \
        template<> \
        struct Deposit<X> {\
            static uint32_t deposit(const std::array<uint32_t, X>& c) { \
                return PImpl::impl_pdep(__VA_ARGS__); \
            } \
        }

ImplSpec(1, c[0]);
ImplSpec(2, c[0], c[1]);
ImplSpec(3, c[0], c[1], c[2]);
ImplSpec(4, c[0], c[1], c[2], c[3]);
ImplSpec(5, c[0], c[1], c[2], c[3], c[4]);

template <unsigned D>
struct Implementation {
    static constexpr unsigned kDimensions = D;

    static uint32_t deposit(const std::array<uint32_t, D>& coords) {
        return Deposit<D>::deposit(coords);
    }

    static std::array<uint32_t, kDimensions> extract(uint32_t cell) {
        return Extract<D>::extract(cell);
    }
    
};

} // namespace BMI2
} // namespace BitManipulationDetails
} // namespace girgs
