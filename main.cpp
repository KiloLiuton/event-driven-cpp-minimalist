#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

#include <omp.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <math.h>
#include <map>
#include "pcg_random/pcg_random.hpp"
#include "dynamics.hpp"
#include "misc.hpp"

#include "experiment-trial_time_evolution.hpp"
#include "experiment-benchmark.hpp"

struct batch_params {
    double coupling_start = 1.0;
    double coupling_end = 3.6;
    int n_batches = 20;
    int batch_id = 0;
    size_t trials = 400;
    size_t iters = 10*N*log(N);
    size_t burn = 10*N*log(N);
    bool verbose = false;
    std::string filename;
};
typedef struct {
    double r;
    double r2;
    double psi;
    double psi2;
    double chi_r;
    double chi_psi;
    float omega;
    double time;
    size_t used_seed;
} Batch;

/* execute a batch of trials and record the average order parameter */
Batch run_batch(
        double coupling,
        size_t trial_iters, size_t trial_burn, size_t trials,
        bool verbose
    );

// PRINTER FUNCTIONS
/* print the state for each site to stdout */
void print_states(States &local_states);
/* print deltas for all sites */
void print_deltas(Deltas &local_deltas);
/* print transitions rates of every site */
void print_rates(Rates &local_rates);
/* get the default file name based on current existing files*/
std::string getDefaultBatchFilename(double a0, double a1, int n);

Batch run_batch(
            double coupling,
            size_t trial_iters, size_t trial_burn, size_t trials,
            bool verbose = false
        ) {
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &current_time);

    // shared variables for each trial
    const size_t seed = 23u * current_time.tv_nsec;
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    Uniform uniform(0.0, 1.0);

    // private variables for each trial
    double r = 0;
    double r2 = 0;
    double psi = 0;
    double psi2 = 0;
    double omega = 0;
    initialize_rates_table(coupling, rates_table);
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
#pragma omp parallel default(none) \
    shared(rates_table) \
    firstprivate(trial_iters,trial_burn,trials,uniform) \
    reduction(+:r,r2,psi,psi2,omega)
    {
        pcg32 RNG(seed);

        States states;
        Deltas deltas;
        Rates rates;
        reset_system(states, deltas, rates, rates_table, RNG);
    #pragma omp for
        for (size_t i = 0; i < trials; i++) {
            pcg32 trial_rng(seed, i);
            initialize_states(states, RNG);
            initialize_deltas(states, deltas);
            initialize_rates(deltas, rates, rates_table);
            Trial trial = run_trial(
                    trial_iters, trial_burn,
                    states, deltas,
                    rates, rates_table,
                    trial_rng,
                    uniform
                );
            r += trial.r;
            r2 += std::pow(trial.r, 2.0);
            psi += trial.psi;
            psi2 += std::pow(trial.psi, 2.0);
            omega += trial.omega;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    Batch batch;
    batch.r = r / trials;
    batch.r2 = r2 / trials;
    batch.psi = psi / trials;
    batch.psi2 = psi2 / trials;
    batch.chi_r = batch.r2 - batch.r * batch.r;
    batch.chi_psi = batch.psi2 - batch.psi * batch.psi;
    batch.omega = omega / trials;
    batch.time = elapsed;
    batch.used_seed = seed;
    return batch;
}

void print_states(States &local_states) {
    for (uint16_t i = 0; i < N; i++) {
        std::cout << unsigned(local_states.array[i]) << ' ';
    }
    std::cout << "(" << local_states.pop[0] << ") "
              << "(" << local_states.pop[1] << ") "
              << "(" << local_states.pop[2] << ")\n";
}

void print_deltas(Deltas &local_deltas) {
    for (uint16_t i = 0; i < N; i++) {
        if (local_deltas[i] >= 0) {
            std::cout << '+' << local_deltas[i] << ' ';
        } else {
            std::cout << local_deltas[i] << ' ';
        }
    }
    int s = 0;
    for (int i = 0; i < N; i++) {
        s += local_deltas[i];
    }
    std::cout << " (sum=" << s << ")\n";
}

void print_rates(Rates &local_rates) {
    std::cout << std::setprecision(2) << std::fixed;
    for (uint16_t i = 0; i < N; i++) {
        std::cout << local_rates.array[i] << ' ';
    }
    std::cout << '\n';
}

std::string getDefaultBatchFilename(double a0, double a1, int n) {
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

int main(int argc, char** argv) {
    std::string help_message =
        "Usage: ./sim[...] [options] [args]\n"
        "-options (args)\n"
        "-h           print this help message\n"
        "--benchmark  run a benchmark for 100 trials\n\n"
        "RUN TRIAL\n"
        "-t           perform a trial with default or specified parameters\n"
        "trial parameters (only available if -t is passed):\n"
        "-tc (real)   [default 2.0] coupling strength for trial\n"
        "-ts (uint)   [default " + std::to_string(DEFAULT_SEED) + "] trial seed\n"
        "-tr (uint)   [default " + std::to_string(DEFAULT_STREAM) + "] trial stream\n"
        "-ti (int)    [default 10*N*log(N)] number of iterations in trial\n"
        "-tb (int)    [default 10*N*log(N)] number of burn iterations in trial\n"
        "-tf (string) trial file name\n\n"
        "RUN BATCH\n"
        "-b           if this flag is present, perform a batch simulation\n"
        "-bs (real)   [default 1.0] initial coupling strength\n"
        "-be (real)   [default 3.6] final coupling strength\n"
        "-bn (int)    [default 20] number of points from initial to final coupling\n"
        "-bt (int)    [default 400] number of trials per coupling value\n"
        "-bi (int)    [default 10*N*log(N)] number of iterations (after burn-in) per\n"
        "             coupling\n"
        "-bb (int)    [default 10*N*log(N)] number of burn-in iterations per coupling\n"
        "-bf (string) batch file name\n";

    // Display help message
    if (
            cmdOptionExists(argv, argv+argc, "-h")
            || cmdOptionExists(argv, argv+argc, "--help")
        ) {
        std::cout << help_message;
        return 0;
    }

    std::cout << "Welcome, to Jurassick Park!\n";
    std::cout << "N=" << N << " K=" << K << " p=" << p
              << " topology_seed=" << TOPOLOGY_SEED << "\n\n";

    // RUN A TRIAL AND LOG IT TO A FILE
    if (cmdOptionExists(argv, argv+argc, "-t")) {
        Time_evolution experiment(argc, argv);
        std::cout << "Logging trial to file: "
                  << experiment.get_filename() << std::endl
                  << "  Coupling: " << experiment.get_coupling() << std::endl
                  << "  Seed: " << experiment.get_seed() << std::endl
                  << "  Stream: " << experiment.get_stream() << std::endl
                  << "  Iters: " << experiment.get_iters() << std::endl
                  << "  Burn: " << experiment.get_burn() << std::endl;
        experiment.run();
    }

    // Get batch parameters
    struct batch_params b_params;
    if (cmdOptionExists(argv, argv+argc, "-b")) {
        try {
            std::string opt;
            if (cmdOptionExists(argv, argv+argc, "-bs")) {
                opt = getCmdOption(argv, argv+argc, "-bs");
                b_params.coupling_start = stof(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-be")) {
                opt = getCmdOption(argv, argv+argc, "-be");
                b_params.coupling_end = stof(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-bn")) {
                opt = getCmdOption(argv, argv+argc, "-bn");
                b_params.n_batches = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-bt")) {
                opt = getCmdOption(argv, argv+argc, "-bt");
                b_params.trials = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-bi")) {
                opt = getCmdOption(argv, argv+argc, "-bi");
                b_params.iters = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-bb")) {
                opt = getCmdOption(argv, argv+argc, "-bb");
                b_params.burn = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-bf")) {
                opt = getCmdOption(argv, argv+argc, "-bf");
                b_params.filename = opt;
            } else {
                b_params.filename = getDefaultBatchFilename(
                        b_params.coupling_start, b_params.coupling_end,
                        b_params.n_batches
                    );
            }
        } catch(...) {
            std::cerr << "Batch parameters not understood!\n"
                         "See -h for help.\n";
            return 0;
        }
    }

    // RUN A BENCHMARK
    if (cmdOptionExists(argv, argv+argc, "--benchmark")) {
        int runs = 20;
        Benchmark benchmark(runs);
        benchmark.run();
        double fastest = benchmark.fastest();
        double slowest = benchmark.slowest();
        double average = benchmark.average();
        std::cout << std::fixed << std::setprecision(6)
                  << "Benchmark results, best of " << runs << " runs.\n"
                  << "  iters: " << (int) 10*N*log(N) << "\n"
                  << "  Min: " << fastest << " s\n"
                  << "  Max: " << slowest << " s\n"
                  << "  Avg: " << average << " s\n";
    }

    // RUN A BATCH OF TRIALS FOR EACH COUPLING STRENGTH
    if (cmdOptionExists(argv, argv+argc, "-b")) {
        std::cout << "\nBatch params:\n"
                  << "coupling_start: " << b_params.coupling_start << '\n'
                  << "coupling_end:   " << b_params.coupling_end << '\n'
                  << "coupling_n:     " << b_params.n_batches << '\n'
                  << "trials:         " << b_params.trials << '\n'
                  << "iters:          " << b_params.iters << '\n'
                  << "burn:           " << b_params.burn << '\n'
                  << "filename:       " << b_params.filename << '\n';
        if (b_params.n_batches < 1) {
            std::cout << "Argument for -bn must be greater than 1!\n";
            return 0;
        }

        FILE* batches_log_file = std::fopen(b_params.filename.c_str(), "w");
        fprintf(
                batches_log_file,
                "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
                "Dynamics_parameters: trials=%lu iters=%lu burn=%lu\n"
                "coupling,r,r2,psi,psi2,chi_r,chi_psi,omega,processing_time,"
                "used_seed\n",
                N, K, p, TOPOLOGY_SEED,
                b_params.trials, b_params.iters, b_params.burn
            );

        for (int i = 0; i < b_params.n_batches; i++) {
            double a;
            if (b_params.n_batches <= 1) {
                a = b_params.coupling_start;
            } else {
                double step = (b_params.coupling_end - b_params.coupling_start)
                              / (b_params.n_batches - 1);
                a = (double) b_params.coupling_start + i * step;
            }
            Batch batch = run_batch(
                    a,
                    b_params.iters,
                    b_params.burn,
                    b_params.trials,
                    b_params.verbose
                );
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
                      << i + 1 << "\\" << b_params.n_batches
                      << "] N=" << N << " K=" << K << " p=" << p << " a=" << a
                      << " Batch finished in: " << std::setprecision(6)
                      << batch.time << "s -- " << std::ctime(&now);
        }
        std::fclose(batches_log_file);
    }

    return 0;
}
