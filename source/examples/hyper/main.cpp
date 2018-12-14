
#include <iostream>
#include <chrono>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>

#include <girgs/girgs-version.h>
#include <girgs/Generator.h>
#include <girgs/Hyperbolic.h>


using namespace std;
using namespace chrono;
using namespace girgs;


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
        clog << "usage: ./hyper\n"
             << "\t\t[-n anInt]          // number of nodes                          default 10000\n"
             << "\t\t[-alpha aFloat]     // ple is 2alpha+1          range [0.5,1]   default 0.75\n"
             << "\t\t[-T aFloat]         // temperature              range [0,1)     default 0\n"
             << "\t\t[-deg anInt]        // average degree           range [1,n)     default 10\n"
             << "\t\t[-rseed anInt]      // radius seed                              default 12\n"
             << "\t\t[-aseed anInt]      // angle seed                               default 130\n"
             << "\t\t[-sseed anInt]      // sampling seed                            default 1400\n"
             << "\t\t[-file aString]     // file name for output graph               default \"graph\"\n"
             << "\t\t[-edge 0|1]         // write result as edgelist (.txt)          default 0\n"
             << "\t\t[-hyp 0|1]          // write hyperbolic coordinates (.hyp)      default 0\n"
             << "\t\t[-stats 0|1]        // shows difference to actual hyp graph     default 0\n";
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ]) : 10000;
    auto alpha  = !params["alpha"].empty()  ? stod(params["alpha"]) : 0.75;
    auto T      = !params["T"    ].empty()  ? stod(params["T"    ]) : 0;
    auto deg    = !params["deg"  ].empty()  ? stoi(params["deg"  ]) : 10;
    auto rseed  = !params["rseed"].empty()  ? stoi(params["rseed"]) : 12;
    auto aseed  = !params["aseed"].empty()  ? stoi(params["aseed"]) : 130;
    auto sseed  = !params["sseed"].empty()  ? stoi(params["sseed"]) : 1400;
    auto file   = !params["file" ].empty()  ? params["file"] : "graph";
    auto edge   = params["edge" ] == "1";
    auto hyp    = params["hyp"  ] == "1";
    auto stats  = params["stats"] == "1";

    // log params and range checks
    cout << "using:\n";
    logParam(n, "n");
    rangeCheck(alpha, 0.5, 1.0, "alpha");
    rangeCheck(T, 0.0, 1.0, "T", false, true);
    rangeCheck(deg, 1, n-1, "deg");
    logParam(rseed, "rseed");
    logParam(aseed, "aseed");
    logParam(sseed, "sseed");
    logParam(file, "file");
    logParam(edge, "edge");
    logParam(hyp, "hyp");
    logParam(stats, "stats");
    auto R = calculateRadius(n, alpha, T, deg);
    logParam(R, "R");
    cout << endl;


    cout << "sampling radii ...\t\t" << flush;
    auto t1 = high_resolution_clock::now();
    auto radii = sampleRadii(n, alpha, R, rseed);
    auto t2 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t2 - t1).count() << "ms"<< endl;

    cout << "sampling angles ...\t\t" << flush;
    auto angles = sampleAngles(n, aseed);
    auto t3 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t3 - t2).count() << "ms" << endl;

    cout << "converting to girg ...\t\t" << flush;
    auto girg_alpha = T == 0 ? std::numeric_limits<double>::infinity() : 1.0 / T;
    auto girg_weights = vector<double>(n, 1.0);
    auto girg_positions = vector<vector<double>>(n, vector<double>(1));
    for(int i = 0; i < n; ++i) {
        girg_weights[i] = radiusToGirgWeight(radii[i], R);
        girg_positions[i][0] = angleToGirgPosition(angles[i]);
    }
    auto t4 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t4 - t3).count() << "ms" << endl;


    cout << "sampling edges for girg ...\t" << flush;
    girgs::Generator generator;
    generator.setWeights(girg_weights);
    generator.setPositions(girg_positions);
    generator.scaleWeights(deg, 1, girg_alpha);
    generator.generate(girg_alpha, sseed);
    auto t5 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t5 - t4).count() << "ms\tavg deg = " << generator.avg_degree() << endl;

    if (edge) {
        cout << "writing edge list (.txt) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        generator.saveEdgeList(file + ".txt");
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    if (hyp) {
        cout << "writing hyp. coords (.hyp) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        auto f = std::ofstream(file);
        for(int i = 0; i < n; ++i)
            f << radii[i] << ' ' << angles[i] << '\n';
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    // compute some stats
    if(!stats)
        return 0;

    cout << "computing stats ...\t\t" << flush;
    auto t6 = high_resolution_clock::now();

    // init accumulators
    long long missing_edges = 0;
    long long false_edges = 0;
    auto missing_deviation = 0.0;
    auto false_deviation = 0.0;

    // normalize adj matrix that edges point from smaller index to larger
    auto adj_list = vector<vector<int>>(n);
    for(auto& node : generator.graph()) {
        for(auto neigh : node.edges) {
            auto mm = minmax(node.index, neigh->index);
            adj_list[mm.first].push_back(mm.second);
        }
    }

    // check all point pairs
    auto current_neighs = std::vector<char>(n, false);
    for(int i=0; i<n; ++i){
        // get one row in adj matrix
        std::fill(current_neighs.begin(), current_neighs.end(), false);
        for(auto neigh : adj_list[i])
            current_neighs[neigh] = true;
        // check that all edges of i are correct in hyperbolic space
        for(int j=i+1; j<n; ++j){
            auto edge_present = current_neighs[j];
            auto dist = girgs::hyperbolicDistance(radii[i], angles[i], radii[j], angles[j]);
            if(dist < R) {
                if(!edge_present) {
                    missing_deviation += R - dist;
                    missing_edges++;
                }
            } else {
                if(edge_present) {
                    false_deviation += dist - R;
                    false_edges++;
                }
            }
        }
    }

    auto t7 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    cout << '\n';
    cout << "missing edges     " << missing_edges << '\n';
    cout << "missing deviation " << missing_deviation << '\n';
    cout << "average mistake   " << missing_deviation / missing_edges << '\n';

    cout << "false edges       " << false_edges << '\n';
    cout << "false deviation   " << false_deviation << '\n';
    cout << "average mistake   " << false_deviation / false_edges << '\n';

    return 0;
}
