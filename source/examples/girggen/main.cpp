
#include <iostream>
#include <chrono>
#include <map>
#include <string>
#include <cstring>

#include <omp.h>

#include <girgs/girgs-version.h>
#include <girgs/Generator.h>


using namespace std;
using namespace chrono;


map<string, string> parseArgs(int argc, char** argv) {
    map<string, string> params;
    for (int i = 1; i < argc; i++) {
        // Get current and next argument
        if (argv[i][0] != '-')
            continue;
        std::string arg = argv[i] + 1; // +1 to skip the -
        // advance one additional position if next is used
        std::string next = (i + 1 < argc ? argv[i++ + 1] : ""); 
        params[std::move(arg)] = std::move(next);
    }
    return params;
}


template<typename T>
void logParam(T value, string name) {
    cout << "\t" << name << "\t=\t" << value << '\n';
}

template<typename T>
void rangeCheck(T value, T min, T max, string name, bool lex = false, bool hex = false) {
    if (value < min || value > max || (value == min && lex) || (value == max && hex)) {
        cerr << "ERROR: parameter " << name << " = " << value << " is not in range "
            << (lex ? "(" : "[") << min << "," << max << (hex ? ")" : "]") << '\n';
        exit(1);
    }
    logParam(value, name);
}



int main(int argc, char* argv[]) {

    // write help
    if (argc < 2 || 0 == strcmp(argv[1], "--help") || 0 == strcmp(argv[1], "-help")) {
        clog << "usage: ./girggen\n"
            << "\t\t[-n anInt]          // number of nodes                          default 10000\n"
            << "\t\t[-d anInt]          // dimension of geometry    range [1,5]     default 1\n"
            << "\t\t[-ple aFloat]       // power law exponent       range (-3,-2]   default -2.5\n"
            << "\t\t[-alpha aFloat]     // model parameter          range (1,inf]   default infinity\n"
            << "\t\t[-deg anInt]        // average degree           range [1,n)     default 10\n"
            << "\t\t[-wseed anInt]      // weight seed                              default 12\n"
            << "\t\t[-pseed anInt]      // position seed                            default 130\n"
            << "\t\t[-sseed anInt]      // sampling seed                            default 1400\n"
            << "\t\t[-threads anInt]    // number of threads to use                 default 1\n"
            << "\t\t[-file aString]     // file name for output graph               default \"graph\"\n"
            << "\t\t[-dot 0|1]          // write result as dot (.dot)               default 0\n"
            << "\t\t[-edge 0|1]         // write result as edgelist (.txt)          default 1\n";
        return 0;
    }

    // write version
    if(argc > 1 && 0 == strcmp(argv[1], "--version")) {
        cout << GIRGS_NAME_VERSION << '\n'
             << GIRGS_PROJECT_DESCRIPTION << '\n'
             << GIRGS_AUTHOR_ORGANIZATION << '\n'
             << GIRGS_AUTHOR_DOMAIN << " (soon)\n"
             << GIRGS_AUTHOR_MAINTAINER << '\n';
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ]) : 10000;
    auto d      = !params["d"    ].empty()  ? stoi(params["d"    ]) : 1;
    auto ple    = !params["ple"  ].empty()  ? stod(params["ple"  ]) : -2.5;
    auto alpha  = !params["alpha"].empty()  ? stod(params["alpha"]) : std::numeric_limits<double>::infinity();
    auto deg    = !params["deg"  ].empty()  ? stoi(params["deg"  ]) : 10;
    auto wseed  = !params["wseed"].empty()  ? stoi(params["wseed"]) : 12;
    auto pseed  = !params["pseed"].empty()  ? stoi(params["pseed"]) : 130;
    auto sseed  = !params["sseed"].empty()  ? stoi(params["sseed"]) : 1400;
    auto threads= !params["threads"].empty()? stoi(params["threads"]) : 1;
    auto file   = !params["file" ].empty()  ? params["file"] : "graph";
    auto dot    = params["dot" ] == "1";
    auto edge   = params["edge"] != "0";

    // log params and range checks
    cout << "using:\n";
    logParam(n, "n");
    rangeCheck(d, 1, 5, "d");
    rangeCheck(ple, -3.0, -2.0, "ple", false, true);
    rangeCheck(alpha, 1.0, std::numeric_limits<double>::infinity(), "alpha", true);
    rangeCheck(deg, 1, n-1, "deg");
    logParam(wseed, "wseed");
    logParam(pseed, "pseed");
    logParam(sseed, "sseed");
    rangeCheck(threads, 1, omp_get_max_threads(), "threads");
    omp_set_num_threads(threads);
    logParam(file, "file");
    logParam(dot, "dot");
    logParam(edge, "edge");
    cout << "\n";

    auto t1 = high_resolution_clock::now();


    cout << "generating weights ...\t\t" << flush;
    girgs::Generator generator;
    generator.setWeights(n, ple, wseed);
    auto t2 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t2 - t1).count() << "ms" << endl;


    cout << "generating positions ...\t" << flush;
    generator.setPositions(n, d, pseed);
    auto t3 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t3 - t2).count() << "ms" << endl;


    cout << "find weight scaling ...\t\t" << flush;
    auto scaling = generator.scaleWeights(deg, d, alpha);
    auto t4 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t4 - t3).count() << "ms\tscaling = " << scaling << endl;

    cout << "sampling edges ...\t\t" << flush;
    generator.generate(alpha, sseed);
    auto t5 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t5 - t4).count() << "ms\tavg deg = " << generator.avg_degree() << endl;

    if (dot) {
        cout << "writing .dot file ...\t\t" << flush;
        auto t6 = high_resolution_clock::now();
        generator.saveDot(file + ".dot");
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    if (edge) {
        cout << "writing edge list (.txt) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        generator.saveEdgeList(file + ".txt");
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    return 0;
}
