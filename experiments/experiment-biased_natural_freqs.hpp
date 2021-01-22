
/*
 *     Perform simulation of a single trial and log the time evolution of all
 * system properties. All simulations parameters are passed on though the
 * `trial_params` struct as well as the topology specified during compilation.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <dynamics.hpp>
#include <misc.hpp>

class Biased_natural_freqs {
public:
    Biased_natural_freqs(int argc, char** argv);
    /* simulate a trial for ITERS time steps after burning BURN steps and save
     * results to log_file. */
    void log_trial_to_file(
            size_t iters, size_t burn,
            States &states, Deltas &deltas,
            Rates &rates, double rates_table[],
            NaturalFreqs &g,
            pcg32 &RNG, Uniform &uniform,
            FILE* log_file
        );
    double run();
    void print_headers();
    double      get_coupling()             { return _t_params.coupling;       }
    int         get_seed()                 { return _t_params.seed;           }
    int         get_stream()               { return _t_params.stream;         }
    int         get_iters()                { return _t_params.iters;          }
    int         get_burn()                 { return _t_params.burn;           }
    std::string get_ic()                   { return _initial_condition;       }
    std::string get_nfreq_bias()           { return _nfreq_bias;              }
    double      get_nfreq_bias_amplitude() { return _nfreq_bias_amplitude;    }
    double      get_timeinterval()         { return _loginterval;             }
    size_t      get_phasebinsize()         { return _phasebinsize;            }
    std::string get_filename()             { return _filename;                }
    bool        get_discardburn()          { return _discardburn;             }
private:
    std::string getDefaultTrialFilename(double coupling);
    void compress_states_and_write_to_file(FILE* file, uint8_t* states, double t);
    const size_t it = 17*N*log(N);
    const size_t bu = 3*N*log(N);
    struct trial_params _t_params {
        2.0, DEFAULT_SEED, DEFAULT_STREAM,
        it, bu
    };
    size_t _phasebinsize = 1+N/100; // create 100 bins by default
    std::string _initial_condition = "wave3";
    std::string _nfreq_bias = "sinesqr";
    double _nfreq_bias_amplitude = 0.5;
    double _loginterval;
    std::string _filename;
    bool _discardburn = false;
};

inline Biased_natural_freqs::Biased_natural_freqs(int argc, char** argv) {
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
                _phasebinsize = stoi(opt);
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
        if (cmdOptionExists(argv, argv+argc, "--nfreqs-bias")) {
            opt = getCmdOption(argv, argv+argc, "--nfreqs-bias");
            if ( (opt=="sinesqr") || (opt=="linear") ) {
                _nfreq_bias = opt;
            }
        }
        if (cmdOptionExists(argv, argv+argc, "--nfreqs-bias-amplitude")) {
            opt = getCmdOption(argv, argv+argc, "--nfreqs-bias-amplitude");
            _nfreq_bias_amplitude = std::stof(opt);
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
        std::cerr << "Experiment parameters not understood!\n"
                     "See -h for help.\n";
        for (int i=0; i<argc; i++) {
            std::cerr << argv[i] << std::endl;
        }
    }
}

inline void Biased_natural_freqs::log_trial_to_file(
            size_t iters, size_t burn,
            States &states, Deltas &deltas, Rates &rates, double rates_table[],
            NaturalFreqs &g,
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
            transition_site(states, deltas, rates, rates_table, g, RNG, uniform);
        }
    } else {
        total_iters = burn + iters;
    }
    double log_counter = std::numeric_limits<float>::infinity();
    for (size_t i=0; i<total_iters; i++) {
        dt = 1.0 / rates.sum;
        time_elapsed += dt;
        transition_site(states, deltas, rates, rates_table, g, RNG, uniform);
        log_counter += dt;
        if (log_counter >= _loginterval) {
            compress_states_and_write_to_file(log_file, states.array, time_elapsed);
            log_counter = 0;
        }
    }
}

inline double Biased_natural_freqs::run() {
    FILE* trial_log_file = fopen((_filename).c_str(), "w");
    fprintf(
            trial_log_file,
            "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
            "Dynamics_parameters: coupling=%f iters=%lu burn=%lu seed=%lu stream=%lu initial_condition=%s discard_burn=%d nfreqs_bias=%s nfreqs_bias_amplitude=%f\n",
            N, K, p, TOPOLOGY_SEED,
            _t_params.coupling, _t_params.iters, _t_params.burn,
            _t_params.seed, _t_params.stream, _initial_condition.c_str(),
            _discardburn, _nfreq_bias.c_str(), _nfreq_bias_amplitude
        );
    fprintf(trial_log_file, "[phases],time_elapsed\n");
    States local_states;
    Deltas local_deltas;
    Rates local_rates;
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    NaturalFreqs g;
    pcg32 RNG(_t_params.seed, _t_params.stream);
    Uniform uniform(0.0, 1.0);

    //  here we manually initialize the natural frequencies array g, instead of
    //  using the `initialize_natural_frequencies` function in order to write
    //  the biased natural frequencies
    if (_nfreq_bias == "sinesqr") {
        for (int n=0; n<N; n++) {
            g[n] = 1.0 + _nfreq_bias_amplitude * (1.0 + sin(n*PI/(2*N))*sin(n*PI/(2*N)));
        }
    } else if (_nfreq_bias == "linear") {
        for (int n=0; n<N/2; n++ ) { g[n] = 1.0 + _nfreq_bias_amplitude * (1.0 + (double) n/N); }
        for (int n=N/2; n>=0; n--) { g[n] = 1.0 + _nfreq_bias_amplitude * (1.0 - (double) n/N); }
    }
    // we then have to manually initialize all other variables related to the
    // dynamics since we can't use the `initialize_everything` function anymore
    initialize_rates_table(_t_params.coupling, rates_table);
    if ( _initial_condition == "random" ) {
        initialize_random_states(local_states, RNG);
    } else if ( _initial_condition == "uniform" ) {
        initialize_uniform_states(local_states);
    } else if ( _initial_condition.substr(0, 4) == "wave" ) {
        initialize_wave_states(local_states, std::stoi(_initial_condition.substr(4)));
    } else {
        std::cout << "Initial condition " << " not understood. Defaulting to random.\n";
        initialize_random_states(local_states, RNG);
    }
    initialize_deltas(local_states, local_deltas);
    initialize_rates(local_deltas, local_rates, rates_table, g);

    struct timespec start, finish;    // measure code run-times
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    log_trial_to_file(
            _t_params.iters, _t_params.burn,
            local_states, local_deltas,
            local_rates, rates_table,
            g,
            RNG, uniform,
            trial_log_file
        );
    clock_gettime(CLOCK_MONOTONIC, &finish);
    fclose(trial_log_file);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    return elapsed;
}

inline void Biased_natural_freqs::print_headers() {
    std::cout
        << "Logging trial to file: " << _filename << '\n'
        << "  Coupling: " << get_coupling() << '\n'
        << "  Seed: " << get_seed() << '\n'
        << "  Stream: " << get_stream() << '\n'
        << "  Iters: " << get_iters() << '\n'
        << "  Burn: " << get_burn() << '\n'
        << "  Discard burn: " << get_discardburn() << '\n'
        << "  Initial condition: " << get_ic() << '\n'
        << "  Log time interval: " << get_timeinterval() << '\n'
        << "  Natural frequency bias: " << get_nfreq_bias() << '\n'
        << "  Natural frequency bias amplitude: " << get_nfreq_bias_amplitude() << '\n'
        << "  Phases compression bin size: " << get_phasebinsize() << '\n';
}

inline std::string Biased_natural_freqs::getDefaultTrialFilename(double coupling) {
    // generate the default file name from the parameters passed to the experiment
    std::ostringstream prefixstream;
    prefixstream.fill('0');
    prefixstream.precision(6);
    prefixstream << "biasedNFtrial-N-";
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

inline void Biased_natural_freqs::compress_states_and_write_to_file(FILE* file, uint8_t* states, double t)
{
    // compress phase values into bins of size `_phasebinsize` and write them to
    // `file` suffixed by the time `t`
    size_t j;
    if (_phasebinsize > 1) {
        size_t compress_counter, n0, n1, n2;
        compress_counter = 0;
        n0 = n1 = n2 = 0;
        for (j=0; j<N; j++) {
            if (states[j] == 0)      n0++;
            else if (states[j] == 1) n1++;
            else if (states[j] == 2) n2++;
            compress_counter++;
            if (compress_counter == _phasebinsize) {
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
    } else if (_phasebinsize == 1) {
        for (j=0; j<N; j++) {
            fprintf(file, "%d,", states[j]);
        }
        fprintf(file, "%f\n", t);
    }
}
