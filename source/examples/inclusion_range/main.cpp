
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


void unifyEdges(Graph& g) {
    for(auto& edge : g)
        if(edge.first > edge.second)
            swap(edge.first, edge.second);
}

// returns the difference of two graphs G1, G2
// G1,2 must be sorted
// first element contains all edges in G1/G2
// second element contains all edges in G2/G1
pair<int, int> edge_diff(const Graph& g1, const Graph& g2) {
    int g1_not_g2 = 0;
    int g2_not_g1 = 0;

    auto it1 = g1.begin();
    auto it2 = g2.begin();
    while(it1!=g1.end() && it2!=g2.end()){
        if(*it1 == *it2){ // edge is in both lists
            it1++;
            it2++;
        } else if(*it1 < *it2) { // edge only in G1
            g1_not_g2++;
            it1++;
        } else { // edge only in G2
            g2_not_g1++;
            it2++;
        }
    }
    g1_not_g2 += g1.end() - it1;
    g2_not_g1 += g2.end() - it2;

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
             << "\t\t[-n1 anInt]         // range start                              default 4096\n"
             << "\t\t[-n2 anInt]         // range end                                default 2097152\n"
             << "\t\t[-alpha aFloat]     // ple is 2alpha+1          range [0.5,1]   default 0.75\n"
             << "\t\t[-deg anInt]        // average degree           range [1,n)     default 100\n"
             << "\t\t[-rseed anInt]      // radius seed                              default 12\n"
             << "\t\t[-aseed anInt]      // angle seed                               default 130\n"
             << "\t\t[-reps anInt]       // repetitions of experiment                default 50\n"
             << "\t\t[-para anInt]       // number of used threads                   default 4\n";
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n1     = !params["n1"   ].empty() ? stoi(params["n1"   ]) : 4096;
    auto n2     = !params["n2"   ].empty() ? stoi(params["n2"   ]) : 2097152;
    auto alpha  = !params["alpha"].empty() ? stod(params["alpha"]) : 0.75; // hrg alpha
    auto deg    = !params["deg"  ].empty() ? stoi(params["deg"  ]) : 100;
    auto rseed  = !params["rseed"].empty() ? stoi(params["rseed"]) : 12;
    auto aseed  = !params["aseed"].empty() ? stoi(params["aseed"]) : 130;
    auto reps   = !params["reps" ].empty() ? stoi(params["reps" ]) : 50;
    auto para  =  !params["para" ].empty() ? stoi(params["para" ]) : 4;

    // log params and range checks
    clog << "using:\n";
    logParam(n1, "n1");
    logParam(n2, "n2");
    rangeCheck(alpha, 0.5, 1.0, "alpha");
    rangeCheck(deg, 1, n1-1, "deg");
    logParam(rseed, "rseed");
    logParam(aseed, "aseed");
    rangeCheck(reps, 0, 100000, "reps");
    rangeCheck(para, 1, 8, "para");
    clog << endl;

    omp_set_num_threads(para);
    auto T = 0.0; // compare two models in threshold versions
    auto sseed = 0; // sampling seed is irrelevant in threshold model

    cout << "n,hrg_deg,low,high,rep\n";
    for(int rep=0; rep < reps; ++rep) {

        // returns required girg deg bounds to be super/sup graph of HRG
        auto deg_low_high = [&](int n) {
            auto R = hypergirgs::calculateRadiusLikeNetworKit(n, alpha, T, deg); // use NetworKit estimator
            // auto R = hypergirgs::calculateRadius(n, hrg_alpha, T, deg);
            auto radii = hypergirgs::sampleRadii(n, alpha, R, rseed+rep+n, false); // this is still done sequential
            auto angles = hypergirgs::sampleAngles(n, aseed+rep+n, false); // this is still done sequential
            auto hrg_edges = hypergirgs::generateEdges(radii, angles, T, R, sseed); // this is parallel
            auto hrg_deg = hrg_edges.size() * 2.0 / n;
            clog << "HRG deg = " << hrg_deg << "\n";
            unifyEdges(hrg_edges);
            sort(hrg_edges.begin(), hrg_edges.end());

            auto girg_alpha = std::numeric_limits<double>::infinity();
            auto girg_weights = vector<double>(n, 1.0);
            auto girg_positions = vector<vector<double>>(n, vector<double>(1));
            for (int i = 0; i < n; ++i) {
                girg_weights[i] = girgs::radiusToGirgWeight(radii[i], R);
                girg_positions[i][0] = girgs::angleToGirgPosition(angles[i]);
            }

            // graph diff of HRG and GIRG
            auto diff = [&](double desired_degree) {
                auto scaled_weights = girg_weights;
                auto scaling = girgs::scaleWeights(scaled_weights, desired_degree, 1, girg_alpha);
                auto girg_edges = girgs::generateEdges(scaled_weights, girg_positions, girg_alpha, sseed); // this is parallel
                unifyEdges(girg_edges);
                sort(girg_edges.begin(), girg_edges.end());
                auto sub = edge_diff(hrg_edges, girg_edges);
                clog << "GIRG  " << desired_degree
                     << "\tdeg " << girg_edges.size() * 2.0 / n
                     << "\tHRG/GIRG " << sub.first
                     << "\tGIRG/HRG " << sub.second << '\n';
                return sub;
            };

            clog << "SCALE UP girg desired avg until girg is a super-graph of HRG\n";
            double low = hrg_deg;
            double high = 2.0 * hrg_deg;
            while (high - low > 0.5) {
                auto mid = (low + high) / 2;
                auto sub = diff(mid);
                if (sub.first == 0)
                    high = mid;
                else
                    low = mid;
            }
            if (low == hrg_deg) cout << "WARNING: default deg no high diff\n";
            auto bound_high = (low + high) / 2;

            clog << "SCALE DOWN girg desired avg until girg is a sub-graph of HRG\n";
            low = 0.5 * hrg_deg;
            high = hrg_deg * 1.2;
            while (high - low > 0.5) {
                auto mid = (low + high) / 2;
                auto sub = diff(mid);
                if (sub.second == 0)
                    low = mid;
                else
                    high = mid;
            }
            if (high == hrg_deg * 1.2) cout << "WARNING: default deg no low diff\n";
            auto bound_low = (low + high) / 2;

            return make_pair(hrg_deg, make_pair(bound_low, bound_high));
        };


        // do actual work
        for (int n = n1; n <= n2; n <<= 1) {
            clog << "Find bounding constants for n= " << n << '\n';

            auto bound_deg = deg_low_high(n);
            cout << n
                 << ',' << bound_deg.first
                 << ',' << bound_deg.second.first
                 << ',' << bound_deg.second.second
                 << ',' << rep
                 << endl;
        }

    } // end reps

    return 0;
}
