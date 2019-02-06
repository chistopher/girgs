#pragma once

#include <algorithm>
#include <vector>
#include <hypergirgs/Hyperbolic.h>

namespace hypergirgs {

class SeedSeq : public std::seed_seq {
public:
    SeedSeq(unsigned int seed, size_t state_size = default_random_engine::state_size)
        : prng_(seed)
        , seeds_(state_size)
    {}

    template< class RandomIt >
    void generate(RandomIt begin, RandomIt end) {
        std::uniform_int_distribution<std::uint32_t> dist;
        std::generate(seeds_.begin(), seeds_.end(),
            [&] {return dist(prng_);});

        std::seed_seq seed_seq(seeds_.begin(), seeds_.end());
        seed_seq.generate(begin, end);
    }

private:
    std::mt19937_64 prng_;
    std::vector<std::uint32_t> seeds_;
};

}