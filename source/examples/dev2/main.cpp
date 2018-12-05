
#include <iostream>
#include <chrono>

#include <girgs/Generator.h>


int main(int argc, char* argv[]) {
    const auto n = 1000000;
    const auto d = 3;
    const auto ple = -2.5;
    //const auto alpha = 1.1;
    const auto alpha = std::numeric_limits<double>::infinity();
    const auto avgDeg = 10;
    const auto weightSeed = 12;
    const auto positionSeed = 130;
    const auto samplingSeed = 1400;

    girgs::Generator generator;
    generator.setPositions(n, d, positionSeed);
    generator.setWeights(n, ple, weightSeed);
    generator.scaleWeights(avgDeg, d, alpha);

    // measure
    auto start = std::chrono::high_resolution_clock::now();
    generator.generate(alpha, samplingSeed);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}
