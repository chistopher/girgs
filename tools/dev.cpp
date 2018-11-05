
#include <iostream>
#include <chrono>

#include <Generator.h>


using namespace std;

int main(int argc, char* argv[])
{

    auto d = 2;
    auto n = 1000;
    auto c = 0.8;
    //auto alpha = 10;
    auto alpha = std::numeric_limits<double>::infinity();
    auto PLE = -2.5;
    auto seed = 15;
    auto weights = generateWeights(n, PLE, seed);
    auto generator = Generator();

    auto start = chrono::high_resolution_clock::now();
    generator.generateGIRG(d, weights, alpha, c, seed);
    auto end = chrono::high_resolution_clock::now();
    cout << d << ":\t" << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;

    cout << generator.avg_degree() << endl;
    
    return 0;
}
