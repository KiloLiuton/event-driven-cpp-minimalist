#define SIN_PHI1 0.8660254037844387
//#ifndef HEADER
//#define HEADER "20-3-0_0-seed_42.h"
//#endif

#include <time.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <math.h>
#include <map>
#include <pcg_random.hpp>

// GLOBAL VARIABLES
struct trial_params {
    double coupling = 2.0;
    size_t iters = 10*N*N;
    size_t burn = 10*N*N;
    size_t seed;
    size_t stream;
    std::string filename;
};
struct batch_params {
    double coupling_start = 1.0;
    double coupling_end = 3.6;
    int coupling_n = 20;
    int batch_id = 0;
    size_t trials = 400;
    size_t iters = 10*N*N;
    size_t burn = 10*N*N;
    std::string filename;
};
struct states {
    uint8_t array[N];
    uint16_t pop[3];
} states;
struct rates {
    double array[N];
    double sum;
} rates;
typedef struct {
    double r;
    double psi;
    double omega;
    double duration;
} trial;
typedef struct {
    double r;
    double r2;
    double psi;
    float omega = 0;
} batch;
int16_t deltas[N];
double ratesTable[NUM_POSSIBLE_TRANSITIONS];
std::uniform_real_distribution<double> UNIFORM(0,1);
pcg32 RNG(42u, 52u);

// INITIALIZER & GETTER FUNCTIONS
/* generates a new random configuration and update all dependencies */
void initialize_everything(
            double coupling,
            size_t seed,
            size_t stream, bool verbose
        );
/* reset the current lattice without changing the coupling strength */
void reset_system(size_t seed, size_t stream);
/* populate the states vector with a random configuration */
void initialize_states();
/* populates delta vector (states must be populated) */
void initialize_deltas();
/* populates rates table (only necessary if coupling changed) */
void initialize_rates_table(double coupling);
/* get the transition rate for a site by accessing deltas vector and then
 * building the access index of the ratesTable. */
double get_rate_from_table(uint16_t site);
/* populates the rates vector */
void initialize_rates();
/* calculates the square of the order parameter for the current state of
 * the system */
double get_squared_op();
/* calculates the order parameter for the current state of the system */
double get_op();
/* calculates the psi order parameter for the current state of the system */
double get_squared_psi_op();
/* calculates the psi order parameter for the current state of the system */
double get_psi_op();

// DYNAMICS FUNCTIONS
/* updates deltas, rates and states for a SINGLE site undergoing transition */
void update_site(int site_index);
/* update rates and delta values for all NEIGHBORS of a site which has been
 * updated */
// TODO remove this function
// void update_neighbors(uint16_t site_index);
/* select an index that will undergo transition */
uint16_t transitionIndex();
/* performs a complete step of the event driven simulation */
uint16_t transition_site();
/* run a trial for ITERS time steps after burning BURN steps and save results
 * to log_file. Logged columns are:
 * r**2,N0,N1,time */
void log_trial_to_file(struct trial_params t_params, FILE* log_file);
/* run a trial and return the average order parameters after BURN iters */
trial run_trial(
        size_t ITERS,
        size_t BURN,
        size_t seed,
        size_t stream
    );
/* run a trial and return the average order parameters and frequency omega
 * after BURN iters */
trial run_omega_trial(struct trial_params t_params);
/* execute a batch of trials and record the average order parameter */
batch run_batch(
            double coupling,
            size_t trials,
            struct trial_params t_params,
            FILE* log_file
        );

// PRINTER FUNCTIONS
/* print the state for each site to stdout */
void print_states();
/* print deltas for all sites */
void print_deltas();
/* print transitions rates of every site */
void print_rates();
/* print the topology included in dynamic header */
void print_geometry();
/* get the default file name based on current existing files*/
std::string getDefaultTrialFilename(double coupling);
/* get the default file name based on current existing files*/
std::string getDefaultBatchFilename(double coupling_start, double coupling_end);
/* print a header to the trial output file (outputs results of a single run) */
void print_file_header(FILE* file, double coupling, int BURN, int ITERS);
/* print a header to the batches output file (many trials for many couplings) */
void print_batches_file_header(FILE* file, int TRIALS, int BURN, int ITERS);

void initialize_everything(double a, size_t seed, size_t stream, bool verbose=false) {
    RNG.seed(seed, stream);
    if (verbose) {
        std::cout << "\nInitializing: "
                  << "N=" << N
                  << "  K=" << K
                  << "  p=" << p
                  << "  a=" << a << '\n';
    }
    initialize_rates_table(a);
    if (verbose) {
        std::cout << "Rates table initialized!\n";
    }
    initialize_states();
    if (verbose) {
        std::cout << "States initialized!\n";
    }
    initialize_deltas();
    if (verbose) {
        std::cout << "Deltas initialized!\n";
    }
    initialize_rates();
    if (verbose) {
        std::cout << "Rates initialized!\n\n";
    }
}

void reset_system(size_t seed, size_t stream) {
    RNG.seed(seed, stream);
    initialize_states();
    initialize_deltas();
    initialize_rates();
}

void initialize_states() {
    states.pop[0] = 0;
    states.pop[1] = 0;
    states.pop[2] = 0;
    for (uint16_t i = 0; i < N; i++) {
        uint8_t state = RNG(3);
        states.array[i] = state;
        switch (state) {
        case 0:
            states.pop[0]++;
            break;
        case 1:
            states.pop[1]++;
            break;
        case 2:
            states.pop[2]++;
            break;
        default:
            std::cout << "Error: Generated state out of range\n";
            break;
        }
    }
}

void initialize_deltas() {
    for (uint16_t i = 0; i < N; i++) {
        const uint16_t ki = NUMBER_OF_NEIGHBORS[i];
        const uint32_t ind = INDEXES[i];
        const int8_t state = states.array[i];
        const int8_t nextState = (state+1)%3;
        // TODO const int8_t nextState = (state + 1 > 2) ? 0 : state + 1;
        int16_t d = 0;
        for (uint16_t j = 0; j < ki; j++) {
            uint16_t nb_ind = NEIGHBOR_LIST[ind+j];
            int8_t nbState = states.array[nb_ind];
            if (nbState == state) {
                d--;
            } else if (nbState == nextState) {
                d++;
            }
        }
        deltas[i] = d;
    }
}

void initialize_rates_table(double a) {
    uint16_t k = K_MIN;
    int16_t d = -K_MIN;
    for (uint32_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        if (d > k) {
            k++;
            d = -k;
        }
        ratesTable[i] = std::exp(a*d/k);
        d++;
    }
}

double get_rate_from_table(uint16_t i) {
    const uint16_t k = NUMBER_OF_NEIGHBORS[i];
    const int16_t d = deltas[i];
    uint32_t idx = (k - K_MIN) * (k + K_MIN) + k + d;
    return ratesTable[idx];
}

void initialize_rates() {
    rates.sum = .0;
    for (uint16_t i = 0; i < N; i++) {
        const double r = get_rate_from_table(i);
        rates.array[i] = r;
        rates.sum += r;
    }
}

double get_squared_op() {
    uint16_t N0 = states.pop[0];
    uint16_t N1 = states.pop[1];
    uint16_t N2 = states.pop[2];
    return (double) (N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / (N*N);
}

double get_op() {
    uint16_t N0 = states.pop[0];
    uint16_t N1 = states.pop[1];
    uint16_t N2 = states.pop[2];
    return (double) std::sqrt(N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / N;
}

double get_squared_psi_op() {
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < N; i++) {
        uint8_t s = states.array[i];
        switch (s) {
        case 0:
            s1 += rates.array[i];
            break;
        case 1:
            s1 += -0.5 * rates.array[i];
            s2 += SIN_PHI1 * rates.array[i];
            break;
        case 2:
            s1 += -0.5 * rates.array[i];
            s2 += -SIN_PHI1 * rates.array[i];
            break;
        default:
            std::cout << (int) s << " <- Unknown state in get_squared_psi_op!\n";
            break;
        }
    }
    return (std::pow(s1, 2.0) + std::pow(s2, 2.0)) / (N*N);
}

double get_psi_op() {
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < N; i++) {
        uint8_t s = states.array[i];
        switch (s) {
        case 0:
            s1 += rates.array[i];
            break;
        case 1:
            s1 += -0.5 * rates.array[i];
            s2 += SIN_PHI1 * rates.array[i];
            break;
        case 2:
            s1 += -0.5 * rates.array[i];
            s2 += -SIN_PHI1 * rates.array[i];
            break;
        default:
            std::cout << (int) s << " <- Unknown state in get_psi_op!\n";
            break;
        }
    }
    return std::sqrt(std::pow(s1, 2.0) + std::pow(s2, 2.0)) / N;
}

void update_site(uint16_t i) {
    uint8_t state = states.array[i];
    uint8_t nextState = (state + 1) % 3;
    // TODO uint8_t nextState = (state + 1 > 2) ? 0 : state + 1;
    states.array[i] = nextState;
    states.pop[state]--;
    states.pop[nextState]++;
    int16_t siteDelta = 0;
    double rateIncrease = 0;
    double rateDecrease = 0;
    for (uint32_t j = INDEXES[i]; j < INDEXES[i] + NUMBER_OF_NEIGHBORS[i]; j++) {
        uint16_t nbIndex = NEIGHBOR_LIST[j];
        uint8_t nbState = states.array[nbIndex];
        if (nbState == state) {
            deltas[nbIndex] += 2;
        } else if (nbState == nextState) {
            siteDelta -= 1;
            deltas[nbIndex] -= 1;
        } else {
            siteDelta += 1;
            deltas[nbIndex] -= 1;
        }
        double nbRate = get_rate_from_table(nbIndex);
        rateIncrease += nbRate;
        rateDecrease -= rates.array[nbIndex];
        rates.array[nbIndex] = nbRate;
    }
    deltas[i] = siteDelta;
    double siteRate = get_rate_from_table(i);
    rateIncrease += siteRate;
    rateDecrease -= rates.array[i];
    rates.array[i] = siteRate;
    rates.sum += rateIncrease + rateDecrease;
}

// TODO remove this function
// void update_neighbors(uint16_t i) {
//     // here we assume that site 'i' has already undergone transition
//     uint8_t nextState = states.array[i];
//     for (uint32_t j = INDEXES[i]; j < INDEXES[i] + NUMBER_OF_NEIGHBORS[i]; j++) {
//         uint16_t nbIndex = NEIGHBOR_LIST[j];
//         uint16_t nbState = states.array[nbIndex];
//         if (nbState == nextState) {
//             deltas[nbIndex] -= 1;
//         } else if (nbState == (nextState + 1) % 3) {
//             deltas[nbIndex] -= 1;
//         } else {
//             deltas[nbIndex] += 2;
//         }
//         // be sure to update rates only after updating `deltas` because
//         // get_rate_from_table uses it to get the updated value for the rate
//         rates.sum -= rates.array[nbIndex];
//         rates.array[nbIndex] = get_rate_from_table(nbIndex);
//         rates.sum += rates.array[nbIndex];
//     }
// }

uint16_t transitionIndex() {
    double partialRate = 0;
    double g = 0;
    double rn = UNIFORM(RNG);
    double randomRate = rn * rates.sum;
    for (uint16_t id = 0; id < N; ++id) {
        g = rates.array[id];
        partialRate += g;
        if (partialRate > randomRate) {
            return id;
        }
    }
    double actualRate =  0;
    for (uint16_t i = 0; i < N; i++) {
        actualRate += rates.array[i];
    }
    std::cout << "Rates refreshed due to Overflow.\n";
    std::cout << "rates.sum: " << rates.sum << " randomRate :" << randomRate << '\n';
    std::cout << "actualRate: " << actualRate << " rn: " << rn << '\n';
    print_rates();

    return N - 1;
}

uint16_t transition_site() {
    uint16_t i = transitionIndex();
    update_site(i);
    return i;
}

void log_trial_to_file(struct trial_params t_params, FILE* log_file) {
    std::cout << "Logging trial to file with:"
              << "  ITERS="    << t_params.iters
              << "  BURN="     << t_params.iters
              << "  coupling=" << t_params.coupling << "\n";

    reset_system(t_params.seed, t_params.stream);

    size_t total_iters = t_params.iters + t_params.burn;
    size_t PROGRESS_INTERVAL = total_iters/20;
    size_t progress_counter = 1;
    double time_elapsed = 0;
    for (size_t i = 0; i < total_iters; ++i) {
        double dt = 1.0 / rates.sum;
        time_elapsed += dt;
        transition_site();

        double r = get_squared_op();
	    double psi = get_psi_op();
        std::fprintf(
                log_file,
                "%16.16f,%16.16f,%d,%d,%f,%f\n",
                r, psi, states.pop[0], states.pop[1], time_elapsed, dt
            );
        if (progress_counter == PROGRESS_INTERVAL) {
            std::cout << std::fixed << std::setprecision(1)
                      << "["
                      << (float) i/total_iters*100
                      << "%]\n";
            progress_counter = 0;
        }
        progress_counter++;
    }
}

trial run_trial(
            size_t I, size_t B,
            size_t seed, size_t stream
        ) {
    RNG.seed(seed, stream);
    initialize_states();
    initialize_deltas();
    initialize_rates();
    for (size_t i = 0; i < B; i++) {
        transition_site();
    }
    double R = 0;
    double PSI = 0;
    double time_elapsed = 0;
    for (size_t i = 0; i < I; i++) {
        double dt = 1.0 / rates.sum;
        R += get_op() * dt;
        PSI += get_psi_op() * dt;
        time_elapsed += dt;
        transition_site();
    }
    trial Trial;
    Trial.r = R / time_elapsed;
    Trial.psi = PSI / time_elapsed;

    return Trial;
}

bool is_crossing(size_t nprev, size_t n, float t, bool is_on_cooldown) {
    if (((nprev <= t && n > t) || (nprev > t && n <= t)) && !is_on_cooldown) {
        return true;
    } else {
        return false;
    }
}

trial run_omega_trial(struct trial_params t_params) {
    for (size_t i = 0; i < t_params.burn; i++) {
        transition_site();
    }

    float t_start = 0;
    float t_end = 0;
    bool started = false;
    bool is_on_cooldown = false;
    int crossings = 0;
    int counter = 0;
    size_t nprev = states.pop[0];
    const int cooldown = N * 1.8;
    const float thold = N / 3.0;

    double R = 0;
    double PSI = 0;
    double time_elapsed = 0;
    for (size_t i = 0; i < t_params.iters; i++) {
        transition_site();
        double dt = 1.0 / rates.sum;
        R += get_op() * dt;
        PSI += get_psi_op() * dt;
        time_elapsed += dt;

        // measure frequency of oscillations by counting how many times the
        // populations cross the N/3 threshold. After detecting such a
        // crossing we need to wait for a cooldown period to expire due to
        // fluctuations.
        size_t n = states.pop[0];
        if (is_crossing(nprev, n, thold, is_on_cooldown)) {
            if (!started) {
                started = true;
                t_start = time_elapsed;
            }
            t_end = time_elapsed;
            is_on_cooldown = true;
            crossings++;
        }
        if (is_on_cooldown) {
            counter++;
        }
        if (counter > cooldown) {
            is_on_cooldown = false;
            counter = 0;
        }
        nprev = n;
    }

    // gather all results and only return a frequency if it exists
    trial Trial;
    Trial.r = R / time_elapsed;
    Trial.psi = PSI / time_elapsed;
    if (crossings > 1) {
        Trial.omega = (crossings - 1) / 2.0 / (t_end - t_start);
    } else {
        Trial.omega = 0;
    }
    Trial.duration = time_elapsed;

    return Trial;
}

batch run_batch(double coupling, size_t trials, struct trial_params t_params, FILE* log_file) {
    initialize_everything(coupling, t_params.seed, 23u, false);
    size_t PROGRESS_INTERVAL = trials / 10;
    size_t progress_counter = 1;
    double r = 0;
    double r2 = 0;
    double psi = 0;
    double omega = 0;
    for (size_t i = 0; i < trials; i++) {
        t_params.seed = i;
        t_params.stream = 2*i;
        trial Trial = run_omega_trial(t_params);
        r += Trial.r;
        r2 += std::pow(Trial.r, 2.0);
        psi += Trial.psi;
        omega += Trial.omega;

        if (progress_counter == PROGRESS_INTERVAL) {
            std::cout << std::setprecision(1) << std::fixed
                      << "[" << (float) i / trials * 100 << "%]\n";
            progress_counter = 0;
        }
        progress_counter++;
    }
    batch Batch;
    Batch.r = r / trials;
    Batch.r2 = r2 / trials;
    Batch.psi = psi / trials;
    Batch.omega = omega / trials;
    return Batch;
}

void print_states() {
    for (uint16_t i = 0; i < N; i++) {
        std::cout << unsigned(states.array[i]) << ' ';
    }
    std::cout << "(" << states.pop[0] << ") "
              << "(" << states.pop[1] << ") "
              << "(" << states.pop[2] << ")\n";
}

void print_deltas() {
    for (uint16_t i = 0; i < N; i++) {
        if (deltas[i] >= 0) {
            std::cout << '+' << deltas[i] << ' ';
        } else {
            std::cout << deltas[i] << ' ';
        }
    }
    int s = 0;
    for (int i = 0; i < N; i++) {
        s += deltas[i];
    }
    std::cout << " (sum=" << s << ")\n";
}

void print_rates() {
    std::cout << std::setprecision(2) << std::fixed;
    for (uint16_t i = 0; i < N; i++) {
        std::cout << rates.array[i] << ' ';
    }
    std::cout << '\n';
}

void print_geometry() {
    std::cout << "INDEXES:\n";
    for (int i = 0; i < N; i++) {
        std::cout << INDEXES[i] << ' ';
    }
    std::cout << '\n';
    std::cout << "NUMBER_OF_NEIGHBORS\n";
    for (int i = 0; i < N; i++) {
        std::cout << NUMBER_OF_NEIGHBORS[i] << ' ';
    }
    std::cout << '\n';
    std::cout << "NEIGHBOR_LIST\n";
    for (int i = 0; i < N*2*K; i++) {
        std::cout << NEIGHBOR_LIST[i] << ' ';
    }
    std::cout << '\n';
}

std::string getDefaultTrialFilename(double coupling) {
    char fname[50];
    sprintf(fname, "N-%05dK-%04dp-%6.6fa-%6.6f_v0.dat", N, K, p, coupling);
    fname[16] = fname[26] = '_';
    int counter = 1;
    while (std::ifstream(fname)) {
        sprintf(
                fname, "N-%05dK-%04dp-%6.6fa-%6.6f_v%d.dat",
                N, K, p, coupling, counter
            );
        fname[16] = '_';
        fname[26] = '_';
        counter++;
    }
    return fname;
}

std::string getDefaultBatchFilename(double coupling_start, double coupling_end) {
    char fname[90];
    sprintf(fname, "batches-N-%05dK-%04dp-%6.6fa-%3.3f-%3.3f_v0.dat",
            N, K, p, coupling_start, coupling_end);
    fname[24] = '_';
    fname[33]  ='_';
    fname[40] = '_';
    int counter = 1;
    while (std::ifstream(fname)) {
        sprintf(fname, "batches-N-%05dK-%04dp-%6.6fa-%3.3f-%3.3f_v%d.dat",
                N, K, p, coupling_start, coupling_end, counter);
        fname[24] = '_';
        fname[33] = '_';
        fname[40] = '_';
        counter++;
    }
    return fname;
}

void print_file_header(FILE* f, double a, int BURN, int ITERS) {
    fprintf(f, "N=%d,k=%d,p=%.6f,a=%.6f,BURN=%d,ITERS=%d\n", N, K, p, a, BURN, ITERS);
    fprintf(f, "r**2,psi,N0,N1,time,dt\n");
}

void print_batches_file_header(FILE* f, int TRIALS, int BURN, int ITERS) {
    fprintf(f, "N=%d, K=%d, p=%f, TRIALS=%d, BURN=%d, ITERS=%d\n", N, K, p, TRIALS, BURN, ITERS);
    fprintf(f, "<<r>>,<<r>^2>,<<psi>>,<<omega>>,a\n");
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::string getCmdOption(char** begin, char** end, const std::string& option) {
    char** s = std::find(begin, end, option);
    if (s != end && ++s != end) {
        return *s;
    }
    return std::string();
}

std::string replaceAllInstances(
            std::string s,
            const std::string& expr,
            const std::string& by
        ) {
    /*Replace all instances of `expr` in s by `by` and return
     * the resulting string*/
    size_t f = s.find(expr);
    if (f == std::string::npos) {
        return s;
    }
    while (f != std::string::npos) {
        s = s.replace(f, expr.length(), by);
        f = s.find(expr);
    }
    return s;
}

int main(int argc, char** argv) {
    std::string help_message = "Usage: ./sim[...] [options] [args]\n"
        "-options (args)\n"
        "-h           print this help message\n\n"
        "-t           perform a trial with default or specified parameters\n"
        "trial parameters (only available if -t is passed):\n"
        "-tc (real)   [default 2.0] coupling strength for trial\n"
        "-tr (uint)   [default 23u] trial seed\n"
        "-ts (uint)   [default 42u] trial stream\n"
        "-ti (int)    [default 10*N*N] number of iterations in trial\n"
        "-tb (int)    [default 10*N*N] number of burn iterations in trial\n"
        "-tf (string) trial file name\n\n"
        "-b           if this flag is present, perform a batch simulation\n"
        "-bs (real)   [default 1.0] initial coupling strength\n"
        "-be (real)   [default 3.6] final coupling strength\n"
        "-bn (real)   [default 20] number of points from initial to final coupling\n"
        "-bt (int)    [default 400] number of trials per coupling value\n"
        "-bi (int)    [default 10*N*N] number of iterations (after burn-in) per\n"
        "             coupling\n"
        "-bb (int)    [default 10*N*N] number of burn-in iterations per coupling\n"
        "-bf (string) batch file name\n";

    // Display help message
    if (cmdOptionExists(argv, argv+argc, "-h")
        || cmdOptionExists(argv, argv+argc, "--help")) {
        std::cout << help_message;
        return 0;
    }

    // Get trial parameters
    struct trial_params t_params;
    if (cmdOptionExists(argv, argv+argc, "-t")) {
        try {
            std::string opt;
            if (cmdOptionExists(argv, argv+argc, "-tc")) {
                opt = getCmdOption(argv, argv+argc, "-tc");
                t_params.coupling = stof(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-tr")) {
                opt = getCmdOption(argv, argv+argc, "-tr");
                t_params.seed = stof(opt);
            } else {
                t_params.seed = 23u;
            }
            if (cmdOptionExists(argv, argv+argc, "-ts")) {
                opt = getCmdOption(argv, argv+argc, "-ts");
                t_params.stream = stof(opt);
            } else {
                t_params.stream = 42u;
            }
            if (cmdOptionExists(argv, argv+argc, "-ti")) {
                opt = getCmdOption(argv, argv+argc, "-ti");
                t_params.iters = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-tb")) {
                opt = getCmdOption(argv, argv+argc, "-tb");
                t_params.burn = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-tf")) {
                opt = getCmdOption(argv, argv+argc, "-tf");
                t_params.filename = opt;
            } else {
                t_params.filename = getDefaultTrialFilename(t_params.coupling);
            }
        } catch(...) {
            std::cerr << "Trial parameters not understood!\n"
                         "See -h for help.\n";
            return 0;
        }
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
                b_params.coupling_n = stoi(opt);
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
                b_params.filename = getDefaultBatchFilename(b_params.coupling_start, b_params.coupling_end);
            }
        } catch(...) {
            std::cerr << "Batch parameters not understood!\n"
                         "See -h for help.\n";
            return 0;
        }
    }

    // variables for measuring code run-times
    struct timespec start, finish;
    double elapsed;

    std::cout << "Welcome, to Jurassick Park!\n";
    std::cout << "N=" << N << " K=" << K << " p=" << p << "\n\n";

    // RUN A TRIAL AND LOG IT TO A FILE
    if (cmdOptionExists(argv, argv+argc, "-t")) {
        std::cout << "\nTrial params:\n";
        std::cout << "coupling: " << t_params.coupling << '\n';
        std::cout << "iters:    " << t_params.iters << '\n';
        std::cout << "burn:     " << t_params.burn << '\n';
        std::cout << "filename: " << t_params.filename << '\n';

        initialize_everything(t_params.coupling, t_params.seed, t_params.stream);
        FILE* trial_log_file = std::fopen(t_params.filename.c_str(), "w");
        clock_gettime(CLOCK_MONOTONIC, &start);
        log_trial_to_file(t_params, trial_log_file);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        std::fclose(trial_log_file);

        // report elapsed time
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        std::cout << std::fixed << std::setprecision(6) << "Elapsed: " << elapsed << "\n";
    }

    // RUN A BATCH OF TRIALS FOR EACH COUPLING STRENGTH
    if (cmdOptionExists(argv, argv+argc, "-b")) {
        std::cout << "\nBatch params:\n";
        std::cout << "coupling start: " << b_params.coupling_start << '\n';
        std::cout << "coupling end:   " << b_params.coupling_end << '\n';
        std::cout << "coupling n:     " << b_params.coupling_n << '\n';
        std::cout << "trials:         " << b_params.trials << '\n';
        std::cout << "iters:          " << b_params.iters << '\n';
        std::cout << "burn:           " << b_params.burn << '\n';
        std::cout << "filename:       " << b_params.filename << '\n';
    }

    return 0;
}
