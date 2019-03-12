#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include "pcg_random/pcg_random.hpp"

struct edge {
    int x, y;
};

class graph {
public:
    typedef std::vector<std::set<int>> kernels;
    graph() {};
    graph(unsigned int n);
    void addnode();
    void addnodes(unsigned int n);
    bool rmnode(int i);
    bool addedge(int i, int j);
    bool rmedge(int i, int j);
    int numspots(unsigned int i);
    bool hasorfan();
    int size() { return _kernels.size(); };
    int numedges() { return _edges; };
    kernels getkernels() { return _kernels; };
    void print();
    void printkernels();
private:
    std::vector<std::set<int>> _kernels;
    int _edges = 0;
};

graph::graph(unsigned int n) {
    addnodes(n);
}

void graph::addnode() {
    _kernels.push_back(std::set<int>());
}

void graph::addnodes(unsigned int n) {
    for (unsigned int i=0; i<n; i++) {
        addnode();
    }
}

bool graph::rmnode(int i) {
    if (i >= size() || i < 0) {
        return false;
    }
    for (int j : _kernels[i]) {
        _kernels[j].erase(i);
        _edges--;
    }
    _kernels.erase(_kernels.begin() + i);
    return true;
}

bool graph::addedge(int i, int j) {
    if (i == j) {
        return false;
    }
    bool p1 = _kernels[i].insert(j).second;
    bool p2 = _kernels[j].insert(i).second;
    _edges++;
    return (p1 && p2);
}

bool graph::rmedge(int i, int j) {
    _edges--;
    return(_kernels[i].erase(j) && _kernels[j].erase(i));
}

int graph::numspots(unsigned int i) {
    return (size() - 1 - _kernels[i].size());
}

bool graph::hasorfan() {
    for (int i=0; i<size(); i++) {
        if (_kernels[i].size() == 0) {
            return false;
        }
    }
    return true;
}

void graph::print() {
    std::set<std::pair<int, int>> seen;
    for (int i=0; i<size(); i++) {
        for (int j : _kernels[i]) {
            int k, l;
            if (j < i) {
                k = j;
                l = i;
            } else {
                k = i;
                l = j;
            }
            if (seen.find(std::pair<int, int> (k, l)) == seen.end()) {
                seen.insert(std::pair<int, int> (k, l));
                std::cout << "(" << k << ", " << l << ")\n";
            }
        }
    }
}

void graph::printkernels() {
    for (int i=0; i<size(); i++) {
        unsigned int s = 0;
        for (int k : _kernels[i]) {
            std::cout << k << ", ";
            s++;
            if (s == _kernels[i].size()) {
                std::cout << std::endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    unsigned int N = std::stoi(argv[1]);
    unsigned int K = std::stoi(argv[2]);
    double p = std::stof(argv[3]);
    unsigned int s;
    if (argc == 5) {
        s = std::stoi(argv[4]);
    } else {
        s = 42u;
    }

    // generate header file name
    std::stringstream sstream;
    sstream << std::fixed << std::setfill('0')
         << std::setw(5) << N << "-"
         << std::setw(4) << K << "-"
         << std::setprecision(6) << p << "-gseed_" << s << ".hpp";
    std::string fname = sstream.str();
    fname = fname.replace(fname.find("."), 1, "_");
    if (std::ifstream(fname)) {
        std::cout << "Header already exists, skipping :)\n";
        return 0;
    }
    std::ofstream header(fname, std::ofstream::out);

    std::uniform_real_distribution<double> uniform(0., 1.);
    pcg32 rng(s);
    graph G = {N};
    for (unsigned int i=1; i<=K; i++) {
        for (unsigned int j=0; j<N; j++) {
            double rn = uniform(rng);
            if (rn < p) {
                if ((G.numspots(j) > 0)) {
                    bool added = false;
                    while (!added) {
                        added = G.addedge(j, rng(N));
                    }
                }
            } else {
                if (j != (j+i)%N) {
                    G.addedge(j, (j+i)%N);
                }
            }
        }
    }

    typedef std::vector<std::set<int>> kernels;
    kernels kers = G.getkernels();
    unsigned int maxk = 0;
    unsigned int mink = N + 1;
    for (std::set<int> k : kers) {
        if (k.size() < mink) mink = k.size();
        if (k.size() > maxk) maxk = k.size();
    }
    unsigned int maxt = (maxk - mink + 1) * (maxk + mink + 1);
    header << "#ifndef TOPOLOGY_HPP\n"
           << "#define TOPOLOGY_HPP\n\n#include <iostream>\n\n"
           << "constexpr uint16_t N = " << N << ";\n"
           << "constexpr uint16_t K = " << K << ";\n"
           << "constexpr float p = " << std::fixed << std::setprecision(6)
               << p << ";\n"
           << "constexpr uint16_t K_MAX = " << maxk << ";\n"
           << "constexpr uint16_t K_MIN = " << mink << ";\n"
           << "constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = " << maxt << ";\n"
           << "constexpr uint32_t TOPOLOGY_SEED = " << s << ";\n"
           << "constexpr uint32_t INDEXES[] = {\n";
    unsigned int i = 0, idx = 0;
    header << idx << ",";
    for (std::set<int> k : kers) {
        i++;
        idx += k.size();
        header << idx << ((i == N-1) ? "\n" : ",");
        if (i == N-1) break;
    }
    header << "};\n"
           << "constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {\n";
    i = 0;
    for (std::set<int> k : kers) {
        i++;
        header << k.size() << ((i == N) ? "\n" : ",");
    }
    header << "};\n"
           << "constexpr uint16_t NEIGHBOR_LIST[] = {\n";
    for (i=0; i<kers.size(); i++) {
        unsigned int x = 1;
        for (int j : kers[i]) {
            header << j << (((i+1 == kers.size()) && (x == kers[i].size())) ? "\n" : ",");
            x++;
        }
    }
    header << "};\n#endif\n";
    // PRINT HEADER

    header.close();
    std::cout << fname << std::endl;
    return 0;
}
