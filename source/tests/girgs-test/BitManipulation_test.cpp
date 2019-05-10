#include <algorithm>
#include <array>
#include <random>
#include <sstream>

// must be here to make ADL work on clang
template<typename T, size_t D>
std::ostream &operator<<(std::ostream &o, const std::array<T, D> &c) {
    std::stringstream ss;
    ss << "[" << c[0];
    for (int d = 1; d < D; d++)
        ss << ", " << c[d];
    ss << "]";
    o << ss.str();
    return o;
}

#include <gtest/gtest.h>
#include <girgs/BitManipulation.h>

template <typename T>
class BitManipulationTest : public ::testing::Test {
public:
    using Implementation = T;
};

using Implementations = ::testing::Types<
#ifdef USE_BMI2
    girgs::BitManipulationDetails::BMI2::Implementation<1>,
    girgs::BitManipulationDetails::BMI2::Implementation<2>,
    girgs::BitManipulationDetails::BMI2::Implementation<3>,
    girgs::BitManipulationDetails::BMI2::Implementation<4>,
    girgs::BitManipulationDetails::BMI2::Implementation<5>,
#endif
    girgs::BitManipulationDetails::Generic::Implementation<1>,
    girgs::BitManipulationDetails::Generic::Implementation<2>,
    girgs::BitManipulationDetails::Generic::Implementation<3>,
    girgs::BitManipulationDetails::Generic::Implementation<4>,
    girgs::BitManipulationDetails::Generic::Implementation<5>
>;

TYPED_TEST_SUITE(BitManipulationTest, Implementations,);

#define COMMON_DEFS \
    using Impl = typename TestFixture::Implementation; \
    constexpr auto D = Impl::kDimensions;              \
    using Coordinate = std::array<uint32_t, D>;        \
    Coordinate coord;                                  \
    std::fill_n(coord.begin(), D, 0);                  \
    std::default_random_engine prng(1 + D);

template <unsigned D>
static uint32_t ReferenceDeposite(std::array<uint32_t, D> coord) {
    uint32_t res = 0;

    for(int i = 0; i < 32; ++i) {
        const auto d = i % D;
        res |= ((coord[d] >> (i / D)) & 1) << i;
    }

    return res;
}

template <unsigned D>
static std::array<uint32_t, D> ReferenceExtract(uint32_t x) {
    std::array<uint32_t, D> res;
    std::fill_n(res.begin(), D, 0);

    for(int i = 0; i < 32; ++i) {
        const auto d = i % D;
        const auto k = i / D;
        res[d] |= ((x >> i) & 1) << k;
    }

    return res;
}



TYPED_TEST(BitManipulationTest, DepositeSingle) {
    COMMON_DEFS

    const auto bits = 32 / D;
    // use 64 bit (at 1ull) to prevent shifts larger or equal to size of type
    const auto max_coord = static_cast<uint32_t>((1ull << bits) - 1);
    std::uniform_int_distribution<unsigned> distr(0, max_coord);

    for(int d = 0; d < D; d++) {

        const auto allowedBits = girgs::BitPattern<D>::kEveryDthBit << d;

        for(int i=0; i < 1000; i++) {
            coord[d] = distr(prng);
            const auto tested = Impl::deposit(coord);
            const auto ref = ReferenceDeposite<D>(coord);

            ASSERT_EQ(tested, tested & allowedBits) << coord;
            ASSERT_EQ(tested, ref) << coord;
        }

        coord[d] = 0; // reset coord to all 0's again
    }
}

TYPED_TEST(BitManipulationTest, DepositeAll) {
    COMMON_DEFS

    const auto bits = 32 / D;
    // use 64 bit (at 1ull) to prevent shifts larger or equal to size of type
    const auto max_coord = static_cast<uint32_t>((1ull << bits) - 1);
    std::uniform_int_distribution<unsigned> distr(0, max_coord);

    std::uniform_int_distribution<unsigned> dim_distr(0, D - 1);
    for(int i=0; i < 10000; i++) {
        const auto d = dim_distr(prng);

        coord[d] = distr(prng);
        const auto tested = Impl::deposit(coord);
        const auto ref = ReferenceDeposite<D>(coord);

        ASSERT_EQ(tested, ref) << coord;
    }
}

TYPED_TEST(BitManipulationTest, ExtractSingle) {
    COMMON_DEFS

    for(int d = 0; d < D; d++) {
        const auto allowedBits = girgs::BitPattern<D>::kEveryDthBit << d;
        std::uniform_int_distribution<unsigned> distr;

        for(int i=0; i < 1000; i++) {
            const auto cell = distr(prng) & allowedBits;

            const auto tested = Impl::extract(cell);
            const auto ref = ReferenceExtract<D>(cell);

            for(int j = 0; j < D; ++j) {
                if (j == d) continue;
                ASSERT_EQ(tested[j], 0) << cell << " tested: " << tested << " ref: " << ref;
            }

            ASSERT_EQ(tested, ref) << cell;
        }
    }
}

TYPED_TEST(BitManipulationTest, ExtractAll) {
    COMMON_DEFS

    std::uniform_int_distribution<unsigned> distr;
    for(int d = 0; d < D; d++) {
        for(int i=0; i < 1000; i++) {
            const auto cell = distr(prng);

            const auto tested = Impl::extract(cell);
            const auto ref = ReferenceExtract<D>(cell);

            ASSERT_EQ(tested, ref) << cell;
        }
    }
}

TYPED_TEST(BitManipulationTest, XRef) {
    COMMON_DEFS

    auto max_bits = (32/D)*D;
    // use 64 bit (at 1ull) to prevent shifts larger or equal to size of type
    auto max_cell = static_cast<uint32_t>((1ull << max_bits) - 1);
    std::uniform_int_distribution<uint32_t> distr(0, max_cell);
    for(int i=0; i < 10000; i++) {
        const auto cell = distr(prng);

        const auto extr = Impl::extract(cell);
        const auto depo = Impl::deposit(extr);

        ASSERT_EQ(cell, depo) << extr;
    }
}
