
#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>

#include <omp.h>

#include <girgs/Generator.h>
#include <girgs/Hyperbolic.h>

#include <hypergirgs/Generator.h>


using namespace std;

using Graph = vector<pair<int,int>>;


void createReverseEdges(Graph& g) {
    const auto n = g.size();
    g.reserve(2*n);
    for(int i=0; i<n; ++i)
        g.push_back({g[i].second, g[i].first});
}

// returns the difference of two graphs G1, G2
// G1,2 must be sorted
// first element contains all edges in G1/G2
// second element contains all edges in G2/G1
pair<Graph, Graph> edge_diff(const Graph& g1, Graph& g2) {
    Graph g1_not_g2;
    Graph g2_not_g1;

    auto it1 = g1.begin();
    auto it2 = g2.begin();
    while(it1!=g1.end() && it2!=g2.end()){
        if(*it1 == *it2){ // edge is in both lists
            it1++;
            it2++;
        } else if(*it1 < *it2) { // edge only in G1
            g1_not_g2.push_back(*it1);
            it1++;
        } else { // edge only in G2
            g2_not_g1.push_back(*it2);
            it2++;
        }
    }
    copy(it1, g1.end(), back_inserter(g1_not_g2));
    copy(it2, g2.end(), back_inserter(g2_not_g1));

    return {g1_not_g2, g2_not_g1};
}


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
    clog << "\t" << name << "\t=\t" << value << '\n';
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
        clog << "Generates input for a HRG and tries to generate a GIRG with the input mapped as described in paper.\n";
        clog << "usage: ./hyper\n"
             << "\t\t[-n anInt]          // number of nodes                          default 10000\n"
             << "\t\t[-alpha aFloat]     // ple is 2alpha+1          range [0.5,1]   default 0.75\n"
             << "\t\t[-deg aFloat]       // average degree           range [1,n)     default 10\n"
             << "\t\t[-odeg aFloat]      // an offset to desired deg for girg        default 0.0\n"
             << "\t\t[-step aFloat]      // odeg step during search                  default 0.1\n"
             << "\t\t[-rseed anInt]      // radius seed                              default 12\n"
             << "\t\t[-aseed anInt]      // angle seed                               default 130\n";
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ]) : 10000;
    auto hrg_alpha  = !params["alpha"].empty()  ? stod(params["alpha"]) : 0.75;
    auto deg    = !params["deg"  ].empty()  ? stod(params["deg"  ]) : 10.0;
    auto odeg   = !params["odeg"].empty()   ? stod(params["odeg" ]) : 0.0;
    auto step   = !params["step"].empty()   ? stod(params["step" ]) : 0.1;
    auto rseed  = !params["rseed"].empty()  ? stoi(params["rseed"]) : 12;
    auto aseed  = !params["aseed"].empty()  ? stoi(params["aseed"]) : 130;

    // log params and range checks
    clog << "using:\n";
    logParam(n, "n");
    rangeCheck(hrg_alpha, 0.5, 1.0, "alpha");
    rangeCheck(deg, 1.0, n-1.0, "deg");
    logParam(odeg, "odeg");
    logParam(step, "step");
    logParam(rseed, "rseed");
    logParam(aseed, "aseed");
    clog << endl;

    omp_set_num_threads(1);
    auto T = 0.0; // compare two models in threshold versions
    auto sseed = 0; // sampling seed is irrelevant in threshold model

    clog << "generate HRG input and edges ";
    auto R = hypergirgs::calculateRadius(n, hrg_alpha, T, deg);
    auto radii = hypergirgs::sampleRadii(n, hrg_alpha, R, rseed, false);
    auto angles = hypergirgs::sampleAngles(n, aseed, false);
    auto hrg_edges = hypergirgs::generateEdges(radii, angles, T, R, sseed);
    clog << "(deg = " << hrg_edges.size()*2.0 / n << ")\n";
    clog << "generated HRG edges: " << hrg_edges.size() << '\n';
    createReverseEdges(hrg_edges);
    sort(hrg_edges.begin(), hrg_edges.end());

    clog << "converting input\n";
    auto girg_alpha = std::numeric_limits<double>::infinity();
    auto girg_weights = vector<double>(n, 1.0);
    auto girg_positions = vector<vector<double>>(n, vector<double>(1));
    for(int i = 0; i < n; ++i) {
        girg_weights[i] = girgs::radiusToGirgWeight(radii[i], R);
        girg_positions[i][0] = girgs::angleToGirgPosition(angles[i]);
    }

    auto diff = [&](double desired_degree) {
        clog << "generate GIRG edges ";
        auto scaled_weights = girg_weights;
        auto scaling = girgs::scaleWeights(scaled_weights, desired_degree, 1, girg_alpha);
        auto girg_edges = girgs::generateEdges(scaled_weights, girg_positions, girg_alpha, sseed);
        clog << "\t(scaling = " << scaling << " deg = " << girg_edges.size()*2.0 / n << ")\n";
        createReverseEdges(girg_edges);
        sort(girg_edges.begin(), girg_edges.end());
        auto sub = edge_diff(hrg_edges, girg_edges);
        clog << "HRG/GIRG " << sub.first.size()/2 << "\tGIRG/HRG " << sub.second.size()/2 << '\n';
        cout << desired_degree << ',' << sub.first.size()/2 << ',' << sub.second.size()/2 << ',' << scaling << ',' << girg_edges.size()/2 << '\n';
        return sub;
    };

    auto high = odeg;
    auto low = odeg-step;

    cout << "deg,hrg_only,girg_only,scaling,girg_edges\n";
    clog << "scale girg desired avg up until girg is a super-graph of HRG\n";
    do {
        auto sub = diff(deg+high);
        if(sub.first.empty()) break;
        high += step;
    } while(true);

    clog << "scale girg desired avg down until girg is a sub-graph of HRG\n";
    do {
        auto sub = diff(deg+low);
        if(sub.second.empty()) break;
        low -= step;
    } while(true);

    return 0;
}
