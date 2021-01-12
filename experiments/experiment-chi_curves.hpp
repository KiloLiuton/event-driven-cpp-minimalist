#include <iomanip>
#include <sstream>
#include <fstream>
#include <chrono>
#include "dynamics.hpp"
#include "misc.hpp"

class Chi_curves {
public:
    Chi_curves(int argc, char** argv);
    void run();
    double get_coupling_start() { return _b_params.coupling_start; }
    double get_coupling_end()   { return _b_params.coupling_end;   }
    double get_num_batches()    { return _b_params.n_batches;      }
    double get_trials()         { return _b_params.trials;         }
    double get_iters()          { return _b_params.iters;          }
    double get_burn()           { return _b_params.burn;           }
    std::string get_ic()        { return _initial_condition;       }
    std::string get_filename()  { return _filename;                }
private:
    size_t it = 17*N*log(N);
    size_t bu = 3*N*log(N);
    struct batch_params _b_params {
        1.0, 3.6, 20, 0, 400, it, bu
    };
    /* get the default file name based on current existing files*/
    std::string getDefaultBatchFilename(double a0, double a1, int n);
    std::string _filename;
    std::string _initial_condition = "uniform";
};

Chi_curves::Chi_curves(int argc, char** argv) {
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
            _b_params.coupling_start = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-be")) {
            opt = getCmdOption(argv, argv+argc, "-be");
            _b_params.coupling_end = stof(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bn")) {
            opt = getCmdOption(argv, argv+argc, "-bn");
            _b_params.n_batches = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bt")) {
            opt = getCmdOption(argv, argv+argc, "-bt");
            _b_params.trials = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bi")) {
            opt = getCmdOption(argv, argv+argc, "-bi");
            _b_params.iters = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bb")) {
            opt = getCmdOption(argv, argv+argc, "-bb");
            _b_params.burn = stoi(opt);
        }
        if (cmdOptionExists(argv, argv+argc, "-bf")) {
            opt = getCmdOption(argv, argv+argc, "-bf");
            _filename = opt;
        } else {
            _filename = getDefaultBatchFilename(
                    _b_params.coupling_start, _b_params.coupling_end,
                    _b_params.n_batches
                );
        }
    } catch(...) {
        std::cerr << "Batch parameters not understood!\n"
                     "See -h for help.\n";
    }
}

void Chi_curves::run() {
    FILE* batches_log_file = std::fopen((_filename).c_str(), "w");
    fprintf(
            batches_log_file,
            "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
            "Dynamics_parameters: trials=%lu iters=%lu burn=%lu "
            "initial_condition=%s\n"
            "coupling,r,r2,psi,psi2,chi_r,chi_psi,omega,processing_time,"
            "used_seed\n",
            N, K, p, TOPOLOGY_SEED,
            _b_params.trials, _b_params.iters, _b_params.burn,
            _initial_condition.c_str()
        );

    for (int i = 0; i < _b_params.n_batches; i++) {
        double a;
        if (_b_params.n_batches <= 1) {
            a = _b_params.coupling_start;
        } else {
            double step = (_b_params.coupling_end - _b_params.coupling_start)
                          / (_b_params.n_batches - 1);
            a = (double) _b_params.coupling_start + i * step;
        }
        Batch batch = run_batch(
                a,
                _b_params.iters,
                _b_params.burn,
                _b_params.trials,
                _initial_condition,
                true
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
                  << i + 1 << "\\" << _b_params.n_batches
                  << "] N=" << N << " K=" << K << " p=" << p << " a=" << a
                  << " Batch finished in: " << std::setprecision(6)
                  << batch.time << "s -- " << std::ctime(&now);
    }
    std::fclose(batches_log_file);
}

std::string Chi_curves::getDefaultBatchFilename(double a0, double a1, int n) {
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
