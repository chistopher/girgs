#pragma once
#include <array>
#include <cstdint>

#ifdef USE_PDEP
#include <immintrin.h>
#endif    


namespace girgs {
namespace bitinterleaving {

namespace word_parallel {
    // fallback for D=3, D >= 5
    template <size_t D>
    struct Impl {
        static uint32_t interleave(const std::array<uint32_t, D>& coords, unsigned targetLevel) {
            unsigned int result = 0u;
            unsigned int bit = 0;

            for(auto l = 0u; l != targetLevel; l++) {
                for(auto d = 0u; d != D; d++) {
                    result |= ((coords[d] >> l) & 1) << bit++;
                }
            }

            return result;
        }
    };

    template <>
    struct Impl<1> {
        static uint32_t interleave(const std::array<uint32_t, 1>& coords, unsigned targetLevel) {
            return coords[0];
        }
    };


    template <>
    struct Impl<2> {
        static uint32_t interleave(const std::array<uint32_t, 2>& coords, unsigned targetLevel) {
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
    struct Impl<4> {
        static uint32_t interleave(const std::array<uint32_t, 4>& coords, unsigned targetLevel) {
            uint64_t z = 0;
            for(int i=0; i != 4; i++)
                z |= static_cast<uint64_t>(coords[i] & 0xff) << (16 * i);

            z = (z | (z << 4)) & 0x0F0F0F0F0F0F0F0F;
            z = (z | (z << 2)) & 0x3333333333333333;
            z = (z | (z << 1)) & 0x5555555555555555;

            z = z | ((z & 0xffff0000ffff0000) >> 15);
            z = z | (z >> 31);

            return static_cast<uint32_t>(z);
        }
    };
} // namespace word_parallel

#ifdef USE_PDEP
namespace PDEP {
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

    // fallback for D > 5
    template <size_t D>
    struct Impl {
        static uint32_t interleave(const std::array<uint32_t, D>& coords, unsigned targetLevel) {
            unsigned int result = 0u;
            unsigned int bit = 0;

            for(auto l = 0u; l != targetLevel; l++) {
                for(auto d = 0u; d != D; d++) {
                    result |= ((coords[d] >> l) & 1) << bit++;
                }
            }

            return result;
        }
    };

    #define ImplSpec(X, ...) \
        template<> \
        struct Impl<X> {\
            static uint32_t interleave(const std::array<uint32_t, X>& c, unsigned = 0) { \
                return PImpl::impl_pdep(__VA_ARGS__); \
            } \
        }

    ImplSpec(1, c[0]);
    ImplSpec(2, c[0], c[1]);
    ImplSpec(3, c[0], c[1], c[2]);
    ImplSpec(4, c[0], c[1], c[2], c[3]);
    ImplSpec(5, c[0], c[1], c[2], c[3], c[4]);

    #undef ImplSpec
} // namespace PDEP
#endif // USE_PDEP

} // namespace bitinterleaving

#ifdef USE_PDEP
template<size_t D>
using BitInterleaving = bitinterleaving::PDEP::Impl<D>;
inline const char* BitInterleavingImpl() {return "PDEP";}
#else
template<size_t D>
using BitInterleaving = bitinterleaving::word_parallel::Impl<D>;
inline const char* BitInterleavingImpl() {return "Word";}
#endif

} // namespace girgs
