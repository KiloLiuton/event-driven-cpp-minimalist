/*
 *     Perform simulation of a single trial and log the time evolution of all
 * system properties. All simulations parameters are passed on though the
 * `trial_params` struct as well as the topology specified during compilation.
 */
#include <sstream>
#include <fstream>
#include "dynamics.hpp"
#include "misc.hpp"

uint8_t median(size_t n0, size_t n1, size_t n2)
{
    // return the phase value of the majority, selecting randomly on ties.
    if (n0 > n2 && n0 > n1)       return 0;
    if (n1 > n2 && n1 > n0)       return 1;
    if (n2 > n0 && n2 > n1)       return 2;
    if ((n0 == n1) && (n1 == n2)) return rand() % 3;
    if (n0 == n1)                 return rand() % 2;
    if (n1 == n2)                 return rand() % 2 + 1;
    if (n0 == n2)                 return (rand() % 2)? 0 : 2;
    std::cout << "oops, this shouldn't have happened\n";
    return 0;
}

class Time_evolution {
public:
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
    void print_headers();
    double get_coupling()      { return _t_params.coupling; }
    int get_seed()             { return _t_params.seed;     }
    int get_stream()           { return _t_params.stream;   }
    int get_iters()            { return _t_params.iters;    }
    int get_burn()             { return _t_params.burn;     }
    std::string get_ic()       { return _initial_condition; }
    double get_timeinterval()  { return _loginterval;       }
    double get_xinterval()     { return _logphases;         }
    std::string get_filename() { return _filename;          }
    bool get_discardburn()     { return _discardburn;       }
private:
    std::string getDefaultTrialFilename(double coupling);
    void compress_states_and_write_to_file(FILE* file, uint8_t* states, double t);
    const size_t it = 17*N*log(N);
    const size_t bu = 3*N*log(N);
    struct trial_params _t_params {
        2.0, DEFAULT_SEED, DEFAULT_STREAM,
        it, bu
    };
    size_t _logphases = 0;
    std::string _initial_condition = "uniform";
    double _loginterval;
    std::string _filename;
    bool _discardburn = false;
};

Time_evolution::Time_evolution(int argc, char** argv) {
    std::string opt;
    try {
        // coupling has to be first since other default values depend on it
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
        if (cmdOptionExists(argv, argv+argc, "--log-phases")) {
            opt = getCmdOption(argv, argv+argc, "--log-phases");
            if (opt.size()) {
                _logphases = stoi(opt);
            } else {
                _logphases = 1;
            }
        }
        if (cmdOptionExists(argv, argv+argc, "--initial-condition")) {
            opt = getCmdOption(argv, argv+argc, "--initial-condition");
            if (opt == "random") {
                _initial_condition = "random";
            } else if (opt == "uniform") {
                _initial_condition = "uniform";
            } else if (opt.substr(0, 4) == "wave") {
                _initial_condition = opt;
            } else {
                std::cout << "Invalid initial condition. Possible values are"
                    " [random, uniform]\nDefaulting to \"" << _initial_condition
                    << "\".\n";
            }
        }
        if (cmdOptionExists(argv, argv+argc, "--log-interval")) {
            opt = getCmdOption(argv, argv+argc, "--log-interval");
            _loginterval = stof(opt);
        } else {
            // double dtmax = N / std::exp(-_t_params.coupling);
            double maxrate = std::exp(std::abs(_t_params.coupling));
            double dtmin = 1.0 / (N * maxrate);
            _loginterval = dtmin * 10;
        }
        if (cmdOptionExists(argv, argv+argc, "-tf")) {
            opt = getCmdOption(argv, argv+argc, "-tf");
            _filename = opt;
        } else {
            _filename = getDefaultTrialFilename(_t_params.coupling);
        }
        if (cmdOptionExists(argv, argv+argc, "--discard-burn")) {
            opt = getCmdOption(argv, argv+argc, "--discard-burn");
            _discardburn = stoi(opt);
        }
    } catch(...) {
        std::cerr << "Trial parameters not understood!\n"
                     "See -h for help.\n";
    }
}

void Time_evolution::log_trial_to_file(
            size_t iters, size_t burn,
            States &states, Deltas &deltas, Rates &rates, double rates_table[],
            pcg32 &RNG, Uniform &uniform,
            FILE* log_file
        ) {
    double time_elapsed = 0;
    double dt;
    size_t total_iters = iters;
    if (_discardburn) {
        for (size_t i=0; i<burn; i++) {
            dt = 1.0 / rates.sum;
            time_elapsed += dt;
            transition_site(states, deltas, rates, rates_table, RNG, uniform);
        }
    } else {
        total_iters = burn + iters;
    }
    double r, psi;
    double log_counter = std::numeric_limits<float>::infinity();
    if (_logphases) {
        for (size_t i=0; i<total_iters; i++) {
            dt = 1.0 / rates.sum;
            time_elapsed += dt;
            transition_site(states, deltas, rates, rates_table, RNG, uniform);
            log_counter += dt;
            if (log_counter >= _loginterval) {
                compress_states_and_write_to_file(log_file, states.array, time_elapsed);
                log_counter = 0;
            }
        }
    } else {
        for (size_t i=0; i<total_iters; i++) {
            dt = 1.0 / rates.sum;
            time_elapsed += dt;
            transition_site(states, deltas, rates, rates_table, RNG, uniform);
            log_counter += dt;
            if (log_counter >= _loginterval) {
                r = get_squared_op(states);
                psi = get_psi_op(states, rates);
                fprintf(log_file,
                        "%16.16f,%16.16f,%d,%d,%f,%f\n",
                        r, psi, states.pop[0], states.pop[1], time_elapsed, dt);
                log_counter = 0;
            }
        }
    }
}

double Time_evolution::run() {
    FILE* trial_log_file = fopen(_filename.c_str(), "w");
    fprintf(
            trial_log_file,
            "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
            "Dynamics_parameters: coupling=%f iters=%lu burn=%lu"
            " seed=%lu stream=%lu initial_condition=%s discard_burn=%d\n",
            N, K, p, TOPOLOGY_SEED,
            _t_params.coupling, _t_params.iters, _t_params.burn,
            _t_params.seed, _t_params.stream, _initial_condition.c_str(),
            _discardburn
        );
    if (_logphases) {
        fprintf(trial_log_file, "[phases],time_elapsed\n");
    } else {
        fprintf(trial_log_file, "r,psi,pop0,pop1,time_elapsed,dt\n");
    }
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
            RNG, _initial_condition,
            false
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

void Time_evolution::print_headers() {
    std::cout
        << "Logging trial to file: " << get_filename() << '\n'
        << "  Coupling: " << get_coupling() << '\n'
        << "  Seed: " << get_seed() << '\n'
        << "  Stream: " << get_stream() << '\n'
        << "  Iters: " << get_iters() << '\n'
        << "  Burn: " << get_burn() << '\n'
        << "  Discard burn: " << get_discardburn() << '\n'
        << "  Initial condition: " << get_ic() << '\n'
        << "  Log time interval: " << get_timeinterval() << '\n';
    if (get_xinterval()) {
        std::cout
            << "  Log phase space interval: " << get_xinterval() << '\n';
    }
}

std::string Time_evolution::getDefaultTrialFilename(double coupling) {
    std::ostringstream prefixstream;
    prefixstream.fill('0');
    prefixstream.precision(6);
    if (_logphases) {
        prefixstream << "phasetrial-N-";
    } else {
        prefixstream << "trial-N-";
    }
    prefixstream.width(5);
    prefixstream << N << "K-";
    prefixstream.width(5);
    prefixstream << K;
    prefixstream << "p-";
    prefixstream << std::fixed << p << "a-" << coupling;
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

void Time_evolution::compress_states_and_write_to_file(FILE* file, uint8_t* states, double t)
{
    size_t j;
    if (_logphases > 1) {
        size_t compress_counter, n0, n1, n2;
        compress_counter = 0;
        n0 = n1 = n2 = 0;
        for (j=0; j<N; j++) {
            if (states[j] == 0)      n0++;
            else if (states[j] == 1) n1++;
            else if (states[j] == 2) n2++;
            compress_counter++;
            if (compress_counter == _logphases) {
                fprintf(file, "%d,", median(n0, n1, n2));
                n0 = n1 = n2 = 0;
                compress_counter = 0;
            }
        }
        if (compress_counter > 0) {
            fprintf(file, "%d,%f\n", median(n0, n1, n2), t);
        } else {
            fprintf(file, "%f\n", t);
        }
    } else if (_logphases == 1) {
        for (j=0; j<N; j++) {
            fprintf(file, "%d,", states[j]);
        }
        fprintf(file, "%f\n", t);
    }
}
