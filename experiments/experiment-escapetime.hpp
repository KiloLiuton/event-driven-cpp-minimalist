/*
 *   Calculate average escape times for various coupling strengths
 * Starting from a uniform IC in state 0, iterate until N0 < N/3 or max iters is
 * reached. This process is repeated a number of times for each coupling
 * strength in the specified range and averages are logged to ouput file.
 */
#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include "dynamics.hpp"
#include "misc.hpp"
#include <omp.h>

class Escape_time {
public:
    Escape_time(int argc, char** argv);
    /* run a trial for ITERS time steps after burning BURN steps and save
     * results to log_file. */
    double run();
    std::vector<double> runescape(double coupling);
    void print_headers();
    std::string get_filename() { return _filename; }
private:
    std::string getDefaultTrialFilename();

    size_t _maxiters = 5*N*log(N);
    size_t _trials = 400;
    std::string _filename;
    double _couplingstart = 4;
    double _couplingend = 4;
    unsigned int _couplingn = 1;
};

Escape_time::Escape_time(int argc, char** argv) {
    std::string opt;
    try {
        // coupling has to be first since other default values depend on it
        if (cmdOptionExists(argv, argv+argc, "-et")) {
            opt = getCmdOption(argv, argv+argc, "-et");
            _trials = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-ei")) {
            opt = getCmdOption(argv, argv+argc, "-ei");
            _maxiters = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-es")) {
            opt = getCmdOption(argv, argv+argc, "-es");
            _couplingstart = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-ee")) {
            opt = getCmdOption(argv, argv+argc, "-ee");
            _couplingend = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-en")) {
            opt = getCmdOption(argv, argv+argc, "-en");
            _couplingn = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-ef")) {
            opt = getCmdOption(argv, argv+argc, "-ef");
            _filename = opt;
        } else {
            _filename = getDefaultTrialFilename();
        }
    } catch(...) {
        std::cerr << "Escape trials parameters not understood!\n"
                     "See -h for help.\n";
    }
}

double Escape_time::run() {
    struct timespec expstart, expfinish;    // measure code run-time
    clock_gettime(CLOCK_MONOTONIC, &expstart);

    FILE* trial_log_file = fopen(_filename.c_str(), "w");
    fprintf(trial_log_file,
            "Graph_parameters: N=%d K=%d p=%f gseed=%d\n"
            "Dynamics_parameters: trials=%lu maxiters=%lu\n",
            N, K, p, TOPOLOGY_SEED, _trials, _maxiters);
    fprintf(trial_log_file, "coupling,[escape_times]\n");
    for (size_t i=0; i<_couplingn; i++) {
        struct timespec trialstart, trialfinish;    // measure code run-time
        clock_gettime(CLOCK_MONOTONIC, &trialstart);
        double a = _couplingstart + i*(_couplingend-_couplingstart)/_couplingn;
        std::vector<double> Tvec = runescape(a);
        fprintf(trial_log_file, "%.8f", a);
        for (auto t: Tvec) {
            fprintf(trial_log_file, ",%.8f", t);
        }
        fprintf(trial_log_file, "\n");
        clock_gettime(CLOCK_MONOTONIC, &trialfinish);
        double trialtime;
        trialtime = (trialfinish.tv_sec - trialstart.tv_sec);
        trialtime += (trialfinish.tv_nsec - trialstart.tv_nsec) / 1000000000.0;

        printf("a=%.4f done (%.2fs) [%lu/%d]\n", a, trialtime, i+1, _couplingn);
    }
    fclose(trial_log_file);

    clock_gettime(CLOCK_MONOTONIC, &expfinish);
    double elapsed;
    elapsed = (expfinish.tv_sec - expstart.tv_sec);
    elapsed += (expfinish.tv_nsec - expstart.tv_nsec) / 1000000000.0;
    return elapsed;
}

std::vector<double> Escape_time::runescape(double a)
{
    std::vector<double> Tvec(_trials, -1);
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    initialize_rates_table(a, rates_table);
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &current_time);
    size_t seed = current_time.tv_nsec;
#pragma omp parallel default(none) \
        shared(rates_table,Tvec,seed)
    {
        Uniform uniform(0.0, 1.0);
        States local_states;
        Deltas local_deltas;
        Rates local_rates;
        NaturalFreqs g;
#pragma omp for
        for (size_t i=0; i<_trials; i++) {
            initialize_uniform_states(local_states);
            initialize_deltas(local_states, local_deltas);
            initialize_rates(local_deltas, local_rates, rates_table, g);
            double t = 0;
            pcg32 RNG(seed, i*17);
            for (size_t j=0; j<_maxiters; j++) {
                transition_site(
                    local_states, local_deltas, local_rates,
                    rates_table, g, RNG, uniform
                );
                 t += 1.0 / local_rates.sum;
                if (local_states.pop[0] <= N/3) {
                    Tvec[i] = t;
                    break;
                }
            } // a normal exit means the system didn't escape and Tvec[i]=-1
        }
    }
    return Tvec;
}

void Escape_time::print_headers() {
    std::cout << "Logging escape times to: " << get_filename() << '\n';
}

std::string Escape_time::getDefaultTrialFilename() {
    std::ostringstream prefixstream;
    prefixstream.fill('0');
    prefixstream.precision(6);
    prefixstream << "escapetrial-N-";
    prefixstream.width(5);
    prefixstream << N << "K-";
    prefixstream.width(5);
    prefixstream << K;
    prefixstream << "p-";
    prefixstream << std::fixed << p << "a-" << _couplingstart;
    prefixstream << "-" << _couplingend << "-" << _couplingn;
    std::string prefix = prefixstream.str();
    prefix = replaceAll(prefix, ".", "_");
    std::ostringstream version;
    version << "_v0.dat";
    std::string fname = prefix + version.str();
    int i = 1;
    while (std::ifstream(fname)) {
        version.str(std::string());  // clear string stream
        version << "_v" << i << ".dat";
        fname = prefix + version.str();
        i++;
    }
    return fname;
}
