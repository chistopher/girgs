
#include <iostream>
#include <chrono>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>

#include <omp.h>

#include <girgs/girgs-version.h>
#include <hypergirgs/Generator.h>


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
void logParam(T value, const string& name) {
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
        clog << "usage: ./genhrg\n"
            << "\t\t[-n anInt]          // number of nodes                          default 10000\n"
            << "\t\t[-alpha aFloat]     // model parameter          range [0.5,1]   default 0.75\n"
            << "\t\t[-t aFloat]         // temperature parameter    range [0,1)     default 0\n"
            << "\t\t[-deg aFloat]       // average degree           range [1,n)     default 10\n"
            << "\t\t[-rseed anInt]      // radii seed                               default 12\n"
            << "\t\t[-aseed anInt]      // angle seed                               default 130\n"
            << "\t\t[-sseed anInt]      // sampling seed                            default 1400\n"
            << "\t\t[-threads anInt]    // number of threads to use                 default 1\n"
            << "\t\t[-nkr 0|1]          // use NetworKit R estimation               default 0\n"
            << "\t\t[-file aString]     // file name for output (w/o ext)           default \"graph\"\n"
            << "\t\t[-edge 0|1]         // write result as edgelist (.txt)          default 0\n"
            << "\t\t[-coord 0|1]        // write hyp. coordinates (.hyp)            default 0\n";
        return 0;
    }

    // write version
    if(argc > 1 && 0 == strcmp(argv[1], "--version")) {
        cout << "HyperGIRGs command line interface.\n\n"
             << GIRGS_NAME_VERSION << '\n'
             << GIRGS_PROJECT_DESCRIPTION << '\n'
             << GIRGS_AUTHOR_ORGANIZATION << '\n'
             << GIRGS_AUTHOR_DOMAIN << '\n'
             << GIRGS_AUTHOR_MAINTAINER << '\n';
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ]) : 10000;
    auto alpha  = !params["alpha"].empty()  ? stod(params["alpha"]) : 0.75;
    auto T      = !params["t"    ].empty()  ? stod(params["t"    ]) : 0;
    auto deg    = !params["deg"  ].empty()  ? stod(params["deg"  ]) : 10.0;
    auto rseed  = !params["rseed"].empty()  ? stoi(params["rseed"]) : 12;
    auto aseed  = !params["aseed"].empty()  ? stoi(params["aseed"]) : 130;
    auto sseed  = !params["sseed"].empty()  ? stoi(params["sseed"]) : 1400;
    auto threads= !params["threads"].empty()? stoi(params["threads"]) : 1;
    auto nkr    = params["nkr"  ] == "1";
    auto file   = !params["file" ].empty()  ? params["file"] : "graph";
    auto edge   = params["edge" ] == "1";
    auto coord  = params["coord"] == "1";

    // log params and range checks
    cout << "using:\n";
    logParam(n, "n");
    rangeCheck(alpha, 0.5, 1.0, "alpha");
    rangeCheck(T, 0.0, 1.0, "t", false, true);
    rangeCheck(deg, 1.0, n-1.0, "deg");
    logParam(rseed, "rseed");
    logParam(aseed, "aseed");
    logParam(sseed, "sseed");
    rangeCheck(threads, 1, omp_get_max_threads(), "threads");
    omp_set_num_threads(threads);
    logParam(nkr, "nkr");
    logParam(file, "file");
    logParam(edge, "edge");
    logParam(coord, "coord");
    cout << "\n";

    cout << "estimate R ...\t\t" << flush;
    auto R = nkr ?
            hypergirgs::calculateRadiusLikeNetworKit(n, alpha, T, deg) :
            hypergirgs::calculateRadius(n, alpha, T, deg);
    cout << "R = " << R << endl;

    auto t1 = high_resolution_clock::now();
    cout << "generating radii ...\t" << flush;
    auto radii = hypergirgs::sampleRadii(n, alpha, R, rseed, threads > 1);
    auto t2 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t2 - t1).count() << "ms" << endl;

    cout << "generating angles ...\t" << flush;
    auto angles = hypergirgs::sampleAngles(n, aseed, threads > 1);
    auto t3 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t3 - t2).count() << "ms" << endl;

    cout << "sampling edges ...\t" << flush;
    auto edges = hypergirgs::generateEdges(radii, angles, T, R, sseed);
    auto t5 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t5 - t3).count() << "ms\tavg deg = " << edges.size()*2.0/n << endl;

    if (edge) {
        cout << "writing edge list (.txt) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        ofstream f{file+".txt"};
        if(!f.is_open()) throw std::runtime_error{"Error: failed to open file \"" + file + ".txt\""};
        f << n << ' ' << edges.size() << "\n\n";
        for(auto& each : edges)
            f << each.first << ' ' << each.second << '\n';
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    if (coord) {
        cout << "writing coordinates (.hyp) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        ofstream f{file+".hyp"};
        if(!f.is_open()) throw std::runtime_error{"Error: failed to open file \"" + file + ".hyp\""};
        f << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        for(auto i{0u}; i<n; ++i)
            f << radii[i]  << ' ' << angles[i] << '\n';
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    return 0;
}
