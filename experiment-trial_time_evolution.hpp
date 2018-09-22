/*
 *     Perform simulation of a single trial and log the time evolution of all
 * system properties. All simulations parameters are passed on though the
 * `trial_params` struct as well as the topology specified during compilation.
 */
#include <sstream>
#include <fstream>
#include "dynamics.hpp"
#include "misc.hpp"

class Time_evolution {
public:
    Time_evolution(struct trial_params t_params);
    Time_evolution(int argc, char** argv);
    /* run a trial for ITERS time steps after burning BURN steps and save
     * results to log_file. */
    void log_trial_to_file(
            size_t iters, size_t burn,
            States &states, Deltas &deltas,
            Rates &rates, double rates_table[],
            pcg32 &RNG, Uniform &uniform,
            FILE* log_file
        );
    double run();
    double get_coupling()      { return _t_params.coupling; }
    int get_seed()             { return _t_params.seed;     }
    int get_stream()           { return _t_params.stream;   }
    int get_iters()            { return _t_params.iters;    }
    int get_burn()             { return _t_params.burn;     }
    std::string get_filename() { return _filename;          }
private:
    const size_t it = 10*N*log(N);
    const size_t bu = 10*N*log(N);
    struct trial_params _t_params {
        2.0, DEFAULT_SEED, DEFAULT_STREAM,
        it, bu
    };
    std::string getDefaultTrialFilename(double coupling);
    std::string _filename;
};

Time_evolution::Time_evolution(struct trial_params t_params) {
    _t_params = t_params;
}
Time_evolution::Time_evolution(int argc, char** argv) {
    try {
        std::string opt;
        if (cmdOptionExists(argv, argv+argc, "-tc")) {
            opt = getCmdOption(argv, argv+argc, "-tc");
            _t_params.coupling = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-ts")) {
            opt = getCmdOption(argv, argv+argc, "-ts");
            _t_params.seed = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-tr")) {
            opt = getCmdOption(argv, argv+argc, "-tr");
            _t_params.stream = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-ti")) {
            opt = getCmdOption(argv, argv+argc, "-ti");
            _t_params.iters = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-tb")) {
            opt = getCmdOption(argv, argv+argc, "-tb");
            _t_params.burn = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-tf")) {
            opt = getCmdOption(argv, argv+argc, "-tf");
            _filename = opt;
        } else {
            _filename = getDefaultTrialFilename(_t_params.coupling);
        }
    } catch(...) {
        std::cerr << "Trial parameters not understood!\n"
                     "See -h for help.\n";
    }
    std::cout << "new constructor called!\n"
              << _t_params.iters << std::endl
              << _t_params.burn << std::endl;
}

void Time_evolution::log_trial_to_file(
            size_t iters, size_t burn,
            States &states, Deltas &deltas, Rates &rates, double rates_table[],
            pcg32 &RNG, Uniform &uniform,
            FILE* log_file
        ) {
    size_t total_iters = iters + burn;
    double time_elapsed = 0;
    for (size_t i = 0; i < total_iters; ++i) {
        double dt = 1.0 / rates.sum;
        time_elapsed += dt;
        transition_site(states, deltas, rates, rates_table, RNG, uniform);

        double r = get_squared_op(states);
	    double psi = get_psi_op(states, rates);
        fprintf(
                log_file,
                "%16.16f,%16.16f,%d,%d,%f,%f\n",
                r, psi, states.pop[0], states.pop[1], time_elapsed, dt
            );
    }
}

double Time_evolution::run() {
    FILE* trial_log_file = fopen(_filename.c_str(), "w");
    fprintf(
            trial_log_file,
            "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
            "Dynamics_parameters: coupling=%f iters=%lu burn=%lu "
            "seed=%lu stream=%lu\n"
            "r,psi,pop0,pop1,time_elapsed,dt\n",
            N, K, p, TOPOLOGY_SEED,
            _t_params.coupling, _t_params.iters, _t_params.burn,
            _t_params.seed, _t_params.stream
        );
    States local_states;
    Deltas local_deltas;
    Rates local_rates;
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    pcg32 RNG(_t_params.seed, _t_params.stream);
    Uniform uniform(0.0, 1.0);
    initialize_everything(
            _t_params.coupling,
            local_states, local_deltas,
            local_rates, rates_table,
            RNG
        );
    struct timespec start, finish;    // measure code run-times
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    log_trial_to_file(
            _t_params.iters, _t_params.burn,
            local_states, local_deltas,
            local_rates, rates_table,
            RNG, uniform,
            trial_log_file
        );
    clock_gettime(CLOCK_MONOTONIC, &finish);
    fclose(trial_log_file);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    return elapsed;
}

std::string Time_evolution::getDefaultTrialFilename(double coupling) {
    std::ostringstream sstream;
    sstream.fill('0');
    sstream.precision(6);
    sstream << "trial-N-";
    sstream.width(5);
    sstream << N << "K-";
    sstream.width(5);
    sstream << K;
    sstream << "p-";
    sstream << std::fixed << p << "a-" << coupling;
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
