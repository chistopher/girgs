
#include <iostream>
#include <chrono>

#include <Generator.h>
#include <DebugHelper.h>


using namespace std;

int main(int argc, char* argv[])
{
    auto weights = generatePowerLawWeights(1000, 1, 1000, -2.1, 1338);
    auto generator = Generator();

    for(auto d=1u; d<4; ++d) {
        auto start = chrono::high_resolution_clock::now();
        generator.generateGIRG(d, weights, 1, 1, 1337);
        auto end = chrono::high_resolution_clock::now();
        cout << d << ":\t" << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
    }

    return 0;
}
