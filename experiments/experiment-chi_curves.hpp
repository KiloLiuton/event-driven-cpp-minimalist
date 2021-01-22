#include <iomanip>
#include <sstream>
#include <fstream>
#include <chrono>
#include <omp.h>
#include "dynamics.hpp"
#include "misc.hpp"

class Chi_curves {
public:
    Chi_curves(int argc, char** argv);
    void run();
    Batch run_batch(double coupling,
                    std::string nfreq_bias,
                    double nfreq_bias_amplitude,
                    std::string ic,
                    size_t trials,
                    size_t iters,
                    size_t burn,
                    int seed);
    double get_coupling_start()       { return coupling_start;        }
    double get_coupling_end()         { return coupling_end;          }
    double get_num_batches()          { return n_batches;             }
    double get_trials()               { return trials;                }
    double get_iters()                { return iters;                 }
    double get_burn()                 { return burn;                  }
    double get_nfreq_bias_amplitude() { return _nfreq_bias_amplitude; }
    std::string get_nfreq_bias()      { return _nfreq_bias;           }
    std::string get_ic()              { return _initial_condition;    }
    std::string get_filename()        { return _filename;             }
private:
    double coupling_start = 1.0;
    double coupling_end = 3.6;
    int n_batches = 20;
    size_t trials = 200;
    size_t iters = 17*N*log(N);
    size_t burn = 3*N*log(N);
    int seed = 23;
    double _nfreq_bias_amplitude = 0.1;
    std::string _nfreq_bias = "sine2";
    /* get the default file name based on current existing files*/
    std::string getDefaultBatchFilename(double a0, double a1, int n);
    std::string _filename;
    std::string _initial_condition = "uniform";
};

inline Chi_curves::Chi_curves(int argc, char** argv) {
    try {
        std::string opt;
        if (cmdOptionExists(argv, argv+argc, "--initial-condition")) {
            opt = getCmdOption(argv, argv+argc, "--initial-condition");
            if (opt != "random" && opt != "uniform") {
                printf("Invalid initial condition. Possible values are "
                       "[random, uniform]\nDefaulting to \"random\".\n");
            } else {
                _initial_condition = opt;
            }
        }
        if (cmdOptionExists(argv, argv+argc, "-bs")) {
            opt = getCmdOption(argv, argv+argc, "-bs");
            coupling_start = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-be")) {
            opt = getCmdOption(argv, argv+argc, "-be");
            coupling_end = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bn")) {
            opt = getCmdOption(argv, argv+argc, "-bn");
            n_batches = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bt")) {
            opt = getCmdOption(argv, argv+argc, "-bt");
            trials = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bi")) {
            opt = getCmdOption(argv, argv+argc, "-bi");
            iters = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bb")) {
            opt = getCmdOption(argv, argv+argc, "-bb");
            burn = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-br")) {
            opt = getCmdOption(argv, argv+argc, "-br");
            seed = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bf")) {
            opt = getCmdOption(argv, argv+argc, "-bf");
            _filename = opt;
        } else {
            _filename = getDefaultBatchFilename(
                    coupling_start, coupling_end,
                    n_batches
                );
        }
    } catch(...) {
        std::cerr << "Batch parameters not understood!\n"
                     "See -h for help.\n";
    }
}

inline void Chi_curves::run() {
    FILE* batches_log_file = std::fopen((_filename).c_str(), "w");
    fprintf(
            batches_log_file,
            "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
            "Dynamics_parameters: trials=%lu iters=%lu burn=%lu "
            "initial_condition=%s\n"
            "coupling,r,r2,psi,psi2,chi_r,chi_psi,omega,processing_time,used_seed\n",
            N, K, p, TOPOLOGY_SEED,
            trials, iters, burn,
            _initial_condition.c_str()
        );

    struct timespec expstart, expfinish;    // measure code run-time
    clock_gettime(CLOCK_MONOTONIC, &expstart);
    for (int i = 0; i < n_batches; i++) {
        double a;
        if (n_batches <= 1) {
            a = coupling_start;
        } else {
            double step = (coupling_end - coupling_start)
                          / (n_batches - 1);
            a = coupling_start + (double) i * step;
        }
        Batch batch = run_batch(a,
                                _nfreq_bias,
                                _nfreq_bias_amplitude,
                                _initial_condition,
                                trials,
                                iters,
                                burn,
                                seed);
        fprintf(
                batches_log_file, 
                "%16.16f,%16.16f,%16.16f,%16.16f,%16.16f,%16.16f,%16.16f,"
                "%16.16f,%16.16f,%lu\n",
                a,
                batch.r,
                batch.r2,
                batch.psi,
                batch.psi2,
                batch.chi_r,
                batch.chi_psi,
                batch.omega,
                batch.time,
                batch.used_seed
           );

        using std::chrono::system_clock;
        std::time_t now = system_clock::to_time_t(system_clock::now());
        std::cout << std::fixed << "["
                  << i + 1 << "\\" << n_batches
                  << "] N=" << N << " K=" << K << " p=" << p << " a=" << a
                  << " Batch finished in: " << std::setprecision(6)
                  << batch.time << "s -- " << std::ctime(&now);
    }
    clock_gettime(CLOCK_MONOTONIC, &expfinish);
    double elapsed = (expfinish.tv_sec - expstart.tv_sec);
    elapsed += (expfinish.tv_nsec - expstart.tv_nsec) / 1000000000.0;
    std::cout << "Chi Curves finished in " << std::setprecision(2) << elapsed << " seconds\n\n";
    std::fclose(batches_log_file);
}

inline Batch Chi_curves::run_batch(double coupling,
                       std::string nfreq_bias,
                       double nfreq_bias_amplitude,
                       std::string ic,
                       size_t trials,
                       size_t iters,
                       size_t burn,
                       int seed)
{
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    initialize_rates_table(coupling, rates_table);

    pcg32 RNG(seed);
    NaturalFreqs g;
    initialize_custom_natural_frequencies(g, nfreq_bias, nfreq_bias_amplitude, RNG);

    size_t i, j;
    double r=0, r2=0, psi=0, psi2=0;
    struct timespec batchstart, batchfinish;    // measure code run-time
    clock_gettime(CLOCK_MONOTONIC, &batchstart);
    omp_set_num_threads(8);
    #pragma omp parallel for default(none) \
        shared(rates_table,g,ic,trials,iters,burn,seed) \
        private(i,j) \
        reduction(+:r,r2,psi,psi2)
    for (i=0; i<trials; i++) {
        States states;
        Deltas deltas;
        Rates rates;
        pcg32 RNG(seed, i*seed);  // trials share a seed value but select different streams
        Uniform uniform(0.0, 1.0);

        if      ( ic=="random"            ) { initialize_random_states(states, RNG);                   }
        else if ( ic=="uniform"           ) { initialize_uniform_states(states);                       }
        else if ( ic.substr(0, 4)=="wave" ) { initialize_wave_states(states, std::stoi(ic.substr(4))); }
        initialize_deltas(states, deltas);
        initialize_rates(deltas, rates, rates_table, g);
        for (j=0; j<burn; j++) {
            transition_site(states, deltas, rates, rates_table, g, RNG, uniform);
        }
        double time_elapsed=0, R=0, PSI=0;
        double dt;
        for (j=0; j<iters; j++) {
            dt = 1.0 / rates.sum;
            time_elapsed += dt;
            transition_site(states, deltas, rates, rates_table, g, RNG, uniform);
            R += get_op(states);
            PSI += get_psi_op(states, rates, g);
        }
        R = R / time_elapsed;
        PSI = PSI / time_elapsed;
        r += R;
        r2 += R*R;
        psi += PSI;
        psi2 += PSI*PSI;
    }
    clock_gettime(CLOCK_MONOTONIC, &batchfinish);
    double elapsed = (batchfinish.tv_sec - batchstart.tv_sec);
    elapsed += (batchfinish.tv_nsec - batchstart.tv_nsec) / 1000000000.0;
    Batch batch;
    batch.r = r / trials;
    batch.r2 = r2 / trials;
    batch.psi = psi / trials;
    batch.psi2 = psi2 / trials;
    batch.chi_r = batch.r2 - batch.r*batch.r;
    batch.chi_psi = batch.psi2 - batch.psi*batch.psi;
    batch.used_seed = seed;
    batch.time = elapsed;
    return batch;
}

inline std::string Chi_curves::getDefaultBatchFilename(double a0, double a1, int n) {
    std::ostringstream sstream;
    sstream.fill('0');
    sstream.precision(6);
    sstream << "batch-N-";
    sstream.width(5);
    sstream << N << "K-";
    sstream.width(5);
    sstream << K;
    sstream << "p-";
    sstream << std::fixed << p << "a-" << a0 << "-" << a1 << "-" << n;
    std::string prefix = sstream.str();
    prefix = replaceAll(prefix, ".", "_");
    std::ostringstream version;
    version << "_v0.dat";
    std::string fname = prefix + version.str();
    int i = 1;
    while (std::ifstream(fname)) {
        version.str(std::string());
        version << "_v" << i << ".dat";
        fname = prefix + version.str();
        i++;
    }
    return fname;
}
