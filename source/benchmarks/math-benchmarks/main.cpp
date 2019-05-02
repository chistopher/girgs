#include <benchmark/benchmark.h>
#include <math.h>
#include <random>

#define POINT_WITH_ORIGINAL
#include <hypergirgs/Generator.h>
#include <hypergirgs/Point.h>

constexpr double kIntToDbl = 1000.;
volatile double dmy = 0.0;

double inv(double x) { return 1.0 / x;}
double pow3(double x)    { return pow(3.0, x); }
double pow100(double x)  { return pow(100.0, x); }
double pow1000(double x) { return pow(1000.0, x);}

#define UNARY_BENCHMARK(X)                         \
    template<typename T>                           \
    static void BM_## X(benchmark::State& state) { \
      T x = state.range(0) / kIntToDbl + dmy;      \
      const T step = 1e10;                         \
      for(auto _ : state) {                        \
        benchmark::DoNotOptimize(x);               \
        const auto tmp = X(x);                     \
        benchmark::DoNotOptimize(tmp);             \
      }                                            \
    }                                              \


#define UNARY_TRIPLE_BENCHMARK(X)               \
    UNARY_BENCHMARK(X)                          \
    BENCHMARK_TEMPLATE(BM_## X, float)          \
        ->Arg(1)                                \
        ->Arg(  1*kIntToDbl)                    \
        ->Arg( 10*kIntToDbl)                    \
        ->Arg(100*kIntToDbl);                   \
    BENCHMARK_TEMPLATE(BM_## X, double)         \
        ->Arg(1)                                \
        ->Arg(  1*kIntToDbl)                    \
        ->Arg( 10*kIntToDbl)                    \
        ->Arg(100*kIntToDbl);                   \


UNARY_TRIPLE_BENCHMARK(inv)
UNARY_TRIPLE_BENCHMARK(pow3)
UNARY_TRIPLE_BENCHMARK(pow100)
UNARY_TRIPLE_BENCHMARK(pow1000)
UNARY_TRIPLE_BENCHMARK(log)
UNARY_TRIPLE_BENCHMARK(log1p)
UNARY_TRIPLE_BENCHMARK(exp)
UNARY_TRIPLE_BENCHMARK(expm1)
UNARY_TRIPLE_BENCHMARK(sqrt)
UNARY_TRIPLE_BENCHMARK(cosh)
UNARY_TRIPLE_BENCHMARK(sinh)
UNARY_TRIPLE_BENCHMARK(acosh)
UNARY_TRIPLE_BENCHMARK(cos)
UNARY_TRIPLE_BENCHMARK(sin)

//////////////////////////////////////////////////////////////////////////////////////////////////////

using Point = hypergirgs::Point;

struct EdgeProbBase {
    EdgeProbBase(double T, double R, double maxProb) : m_T(T), m_R(R), m_maxProb(maxProb) {}

    bool operator() (const Point& a, const Point& b, double uni_rnd) const noexcept {
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
        benchmark::DoNotOptimize(uni_rnd);
        return false;
    }

    double m_T;
    double m_R;
    double m_maxProb;
};

struct EdgeProbNaive : public EdgeProbBase {
    EdgeProbNaive(double T, double R, double maxProb) : EdgeProbBase(T,R,maxProb) {}
    bool operator() (const Point& a, const Point& b, double uni_rnd) const noexcept {
        const auto dist = acosh(std::max(1.,
            cosh(a.radius - b.radius) + (1. - cos(a.angle - b.angle)) * sinh(a.radius) * sinh(b.radius)));
        const auto connection_prob = 1.0 / (1.0 + std::exp(0.5/m_T*(dist-m_R)));
        return uni_rnd < connection_prob / m_maxProb;
    }
};

struct EdgeProbSimplified : public EdgeProbBase {
    EdgeProbSimplified(double T, double R, double maxProb) : EdgeProbBase(T,R,maxProb), m_scale(0.5 / T) {}
    bool operator() (const Point& a, const Point& b, double uni_rnd) const noexcept {
        const auto dist = acosh(std::max(1., (a.coth_r * b.coth_r - a.cos_phi * b.cos_phi - a.sin_phi * b.sin_phi) / (a.invsinh_r * b.invsinh_r)));
        const auto connection_prob = (1.0 + std::exp(m_scale*(dist-m_R)));
        return uni_rnd * m_maxProb * connection_prob < 1.0;
    }

    double m_scale;
};


template <typename Impl>
static void BM_edge_prob(benchmark::State& state) {
    std::mt19937_64 prng;
    const size_t n = 10000;
    const double R = 20;

    std::uniform_real_distribution<double> dist_prob;

    Impl base{2.0, R, 1e-5};

    std::vector<hypergirgs::Point> points;
    points.reserve(n);
    {
        const auto angles = hypergirgs::sampleAngles(n, 1, false);
        const auto radii = hypergirgs::sampleRadii(n, 1.0, R, false);
        for(int i=0; i != n; i++)
            points.emplace_back(i, radii[i], angles[i]);
    }

    int i = 0;
    for(auto _ : state) {
        i = i < n - 2 ? i + 2 : 0;

        const auto is_edge = base(points[i], points[i+1], dist_prob(prng));
        benchmark::DoNotOptimize(is_edge);
    }
}

BENCHMARK_TEMPLATE(BM_edge_prob, EdgeProbBase);
BENCHMARK_TEMPLATE(BM_edge_prob, EdgeProbNaive);
BENCHMARK_TEMPLATE(BM_edge_prob, EdgeProbSimplified);



BENCHMARK_MAIN();

//         x += step;                                 \