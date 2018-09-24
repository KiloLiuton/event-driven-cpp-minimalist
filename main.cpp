#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

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
#include "pcg_random/pcg_random.hpp"
#include "dynamics.hpp"
#include "misc.hpp"

#include "experiment-trial_time_evolution.hpp"
#include "experiment-benchmark.hpp"
#include "experiment-chi_curves.hpp"

// PRINTER FUNCTIONS
/* print the state for each site to stdout */
void print_states(States &local_states);
/* print deltas for all sites */
void print_deltas(Deltas &local_deltas);
/* print transitions rates of every site */
void print_rates(Rates &local_rates);


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

int main(int argc, char** argv) {
    std::string help_message =
        "Usage: ./sim[...] [options] [args]\n"
        "-options (args)\n"
        "-h           print this help message\n"
        "--benchmark  run a benchmark\n\n"
        "TRIAL TIME EVOLUTION\n"
        "-t           perform a trial with default or specified parameters\n"
        "trial parameters (only available if -t is passed):\n"
        "-tc (real)   [default 2.0] coupling strength for trial\n"
        "-ts (uint)   [default " + std::to_string(DEFAULT_SEED) + "] trial seed\n"
        "-tr (uint)   [default " + std::to_string(DEFAULT_STREAM) + "] trial stream\n"
        "-ti (int)    [default 10*N*log(N)] number of iterations in trial\n"
        "-tb (int)    [default 10*N*log(N)] number of burn iterations in trial\n"
        "-tf (string) trial file name\n\n"
        "CHI CURVE\n"
        "-b           generate a chi curve for the given coupling values\n"
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

    // RUN A BATCH OF TRIALS FOR EACH COUPLING STRENGTH
    if (cmdOptionExists(argv, argv+argc, "-b")) {
        Chi_curves exp(argc, argv);
        std::cout << "Chi curves parameters:\n"
                  << "coupling_start: " << exp.get_coupling_start() << '\n'
                  << "coupling_end:   " << exp.get_coupling_end() << '\n'
                  << "coupling_n:     " << exp.get_num_batches() << '\n'
                  << "trials:         " << exp.get_trials() << '\n'
                  << "iters:          " << exp.get_iters() << '\n'
                  << "burn:           " << exp.get_burn() << '\n'
                  << "filename:       " << exp.get_filename() << '\n';
        if (exp.get_num_batches() < 1) {
            std::cout << "Argument for -bn must be greater than 1!\n";
            return 0;
        }
        exp.run();
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

    return 0;
}
