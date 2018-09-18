#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

#include <omp.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <math.h>
#include <map>
#include <pcg_random.hpp>

struct trial_params {
    double coupling = 2.0;
    size_t seed = DEFAULT_SEED;
    size_t stream = DEFAULT_STREAM;
    size_t iters = 10*N*log(N);
    size_t burn = 10*N*log(N);
    std::string filename;
};
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
    uint8_t array[N];
    uint16_t pop[3];
} States;
typedef int16_t Deltas[N];
typedef struct {
    double array[N];
    double sum;
} Rates;
typedef struct {
    double r;
    double psi;
    double omega;
    double duration;
} Trial;
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

// GLOBAL VARIABLES
typedef std::uniform_real_distribution<double> Uniform;
double ratesTable[NUM_POSSIBLE_TRANSITIONS];
// Deltas deltas;

// INITIALIZER & GETTER FUNCTIONS
/* generates a new random configuration and update all dependencies */
void initialize_everything(
        double coupling,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        bool verbose
    );
/* reset the current lattice without changing the coupling strength */
void reset_system(
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG
    );
/* populates rates table with a given coupling value */
void initialize_rates_table(double coupling);
/* get the transition rate based on the current delta value of site i */
double get_rate_from_table(uint16_t i, int16_t d);
/* populate a given states vector with a random configuration */
void initialize_states(States &local_states, pcg32 &RNG);
/* populates a given rates vector from the current rates table */
void initialize_rates(Deltas &local_deltas, Rates &local_rates);
/* populates delta vector (states must be populated) */
void initialize_deltas(States &local_states, Deltas &local_deltas);
/* calculate the squared order parameter of a given states array */
double get_squared_op(States &local_states);
/* calculates the order parameter of a given states array */
double get_op(States &local_states);
/* calculates the psi order parameter for the current state of the system */
double get_squared_psi_op(States &local_states, Rates &local_rates);
/* calculates the psi order parameter for the current state of the system */
double get_psi_op(States &local_states, Rates &local_rates);

// DYNAMICS FUNCTIONS
/* updates deltas, rates and states for a site and its neighbors */
void update_site(
        int site_index,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates
    );
/* select an index that will undergo transition */
uint16_t transitionIndex(
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    );
/* performs a complete step of the event driven simulation */
uint16_t transition_site(
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    );
/* run a trial for ITERS time steps after burning BURN steps and save results
 * to log_file. Logged columns are:
 * r**2,N0,N1,time */
void log_trial_to_file(
        size_t iters,
        size_t burn,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform,
        FILE* log_file,
        bool verbose
    );
/* run a trial and return the average order parameters after BURN iters */
Trial run_no_omega_trial(
        size_t iters,
        size_t burn,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    );
/* run a trial and return the average order parameters and frequency omega */
Trial run_trial(
        size_t iters,
        size_t burn,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    );
/* determine if there is a crossing of threshold for frequency detection */
bool is_crossing(size_t nprev, size_t n, float t, bool is_on_cooldown);
/* execute a batch of trials and record the average order parameter */
Batch run_batch(
        double coupling,
        size_t trial_iters,
        size_t trial_burn,
        size_t trials,
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
std::string getDefaultTrialFilename(double coupling);
/* get the default file name based on current existing files*/
std::string getDefaultBatchFilename(double coupling_start, double coupling_end);

void initialize_everything(
        double a,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        bool verbose=false
    ) {
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
    initialize_states(local_states, RNG);
    if (verbose) {
        std::cout << "States initialized!\n";
    }
    initialize_deltas(local_states, local_deltas);
    if (verbose) {
        std::cout << "Deltas initialized!\n";
    }
    initialize_rates(local_deltas, local_rates);
    if (verbose) {
        std::cout << "Rates initialized!\n\n";
    }
}

void reset_system(
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG
    ) {
    initialize_states(local_states, RNG);
    initialize_deltas(local_states, local_deltas);
    initialize_rates(local_deltas, local_rates);
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

double get_rate_from_table(uint16_t i, int16_t d) {
    const uint16_t k = NUMBER_OF_NEIGHBORS[i];
    uint32_t idx = (k - K_MIN) * (k + K_MIN) + k + d;
    return ratesTable[idx];
}

void initialize_states(States &local_states, pcg32 &RNG) {
    local_states.pop[0] = 0;
    local_states.pop[1] = 0;
    local_states.pop[2] = 0;
    for (uint16_t i = 0; i < N; i++) {
        uint8_t state = RNG(3);
        local_states.array[i] = state;
        switch (state) {
        case 0:
            local_states.pop[0]++;
            break;
        case 1:
            local_states.pop[1]++;
            break;
        case 2:
            local_states.pop[2]++;
            break;
        default:
            std::cout << "Error: Generated state out of range\n";
            break;
        }
    }
}

void initialize_rates(Deltas &local_deltas, Rates &local_rates) {
    local_rates.sum = .0;
    for (uint16_t i = 0; i < N; i++) {
        int d = local_deltas[i];
        const double r = get_rate_from_table(i, d);
        local_rates.array[i] = r;
        local_rates.sum += r;
    }
}

void initialize_deltas(States &local_states, Deltas &local_deltas) {
    for (uint16_t i = 0; i < N; i++) {
        const uint16_t ki = NUMBER_OF_NEIGHBORS[i];
        const uint32_t ind = INDEXES[i];
        const int8_t state = local_states.array[i];
        // TODO
        // const int8_t nextState = (state+1) % 3;
        const int8_t nextState = (state + 1 > 2) ? 0 : state + 1;
        int16_t d = 0;

        for (uint16_t j = 0; j < ki; j++) {
            uint16_t nb_ind = NEIGHBOR_LIST[ind+j];
            int8_t nbState = local_states.array[nb_ind];
            if (nbState == state) {
                d--;
            } else if (nbState == nextState) {
                d++;
            }
        }
        local_deltas[i] = d;
    }
}

double get_squared_op(States &local_states) {
    uint16_t N0 = local_states.pop[0];
    uint16_t N1 = local_states.pop[1];
    uint16_t N2 = local_states.pop[2];
    return (double) (N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / (N*N);
}

double get_op(States &local_states) {
    uint16_t N0 = local_states.pop[0];
    uint16_t N1 = local_states.pop[1];
    uint16_t N2 = local_states.pop[2];
    return (double) sqrt(N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / N;
}

double get_squared_psi_op(States &local_states, Rates &local_rates) {
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < N; i++) {
        uint8_t s = local_states.array[i];
        switch (s) {
        case 0:
            s1 += local_rates.array[i];
            break;
        case 1:
            s1 += -0.5 * local_rates.array[i];
            s2 += SIN_PHI1 * local_rates.array[i];
            break;
        case 2:
            s1 += -0.5 * local_rates.array[i];
            s2 += -SIN_PHI1 * local_rates.array[i];
            break;
        default:
            std::cout << (int) s << "<- Unknown state in get_squared_psi_op!\n";
            break;
        }
    }
    return (std::pow(s1, 2.0) + std::pow(s2, 2.0)) / (N*N);
}

double get_psi_op(States &local_states, Rates &local_rates) {
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < N; i++) {
        uint8_t s = local_states.array[i];
        switch (s) {
       case 0:
            s1 += local_rates.array[i];
            break;
        case 1:
            s1 += -0.5 * local_rates.array[i];
            s2 += SIN_PHI1 * local_rates.array[i];
            break;
        case 2:
            s1 += -0.5 * local_rates.array[i];
            s2 += -SIN_PHI1 * local_rates.array[i];
            break;
        default:
            std::cout << (int) s << " <- Unknown state in get_psi_op!\n";
            break;
        }
    }
    return std::sqrt(std::pow(s1, 2.0) + std::pow(s2, 2.0)) / N;
}

void update_site(
        uint16_t i,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates
    ) {
    uint8_t state = local_states.array[i];
    // TODO
    // uint8_t nextState = (state + 1) % 3;
    uint8_t nextState = (state + 1 > 2) ? 0 : state + 1;
    local_states.array[i] = nextState;
    local_states.pop[state]--;
    local_states.pop[nextState]++;
    int16_t siteDelta = 0;
    double rateIncrease = 0;
    double rateDecrease = 0;
    uint16_t ki = NUMBER_OF_NEIGHBORS[i];
    for (uint32_t j = INDEXES[i]; j < INDEXES[i] + ki; j++) {
        uint16_t nbIndex = NEIGHBOR_LIST[j];
        uint8_t nbState = local_states.array[nbIndex];
        if (nbState == state) {
            local_deltas[nbIndex] += 2;
        } else if (nbState == nextState) {
            siteDelta -= 1;
            local_deltas[nbIndex] -= 1;
        } else {
            siteDelta += 1;
            local_deltas[nbIndex] -= 1;
        }
        double nbRate = get_rate_from_table(nbIndex, local_deltas[nbIndex]);
        rateIncrease += nbRate;
        rateDecrease -= local_rates.array[nbIndex];
        local_rates.array[nbIndex] = nbRate;
    }
    local_deltas[i] = siteDelta;
    double siteRate = get_rate_from_table(i, local_deltas[i]);
    rateIncrease += siteRate;
    rateDecrease -= local_rates.array[i];
    local_rates.array[i] = siteRate;
    local_rates.sum += rateIncrease + rateDecrease;
}

uint16_t transitionIndex(
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    ) {
    double partialRate = 0;
    double g = 0;
    double rn = uniform(RNG);
    double randomRate = rn * local_rates.sum;
    for (uint16_t id = 0; id < N; ++id) {
        g = local_rates.array[id];
        partialRate += g;
        if (partialRate > randomRate) {
            return id;
        }
    }
    double actualRate =  0;
    for (uint16_t i = 0; i < N; i++) {
        actualRate += local_rates.array[i];
    }
    std::cout << "Rates refreshed due to Overflow.\n";

    return N - 1;
}

uint16_t transition_site(
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    ) {
    uint16_t i = transitionIndex(local_rates, RNG, uniform);
    update_site(i, local_states, local_deltas, local_rates);
    return i;
}

void log_trial_to_file(
            size_t iters,
            size_t burn,
            States &local_states,
            Deltas &local_deltas,
            Rates &local_rates,
            pcg32 &RNG,
            Uniform &uniform,
            FILE* log_file
        ) {
    size_t total_iters = iters + burn;
    double time_elapsed = 0;
    for (size_t i = 0; i < total_iters; ++i) {
        double dt = 1.0 / local_rates.sum;
        time_elapsed += dt;
        transition_site(local_states, local_deltas, local_rates, RNG, uniform);

        double r = get_squared_op(local_states);
	    double psi = get_psi_op(local_states, local_rates);
        fprintf(
                log_file,
                "%16.16f,%16.16f,%d,%d,%f,%f\n",
                r, psi, local_states.pop[0], local_states.pop[1],
                time_elapsed, dt
            );
    }
}

Trial run_no_omega_trial(
        size_t iters,
        size_t burn,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    ) {
    initialize_states(local_states, RNG);
    initialize_deltas(local_states, local_deltas);
    initialize_rates(local_deltas, local_rates);
    for (size_t i = 0; i < burn; i++) {
        transition_site(local_states, local_deltas, local_rates, RNG, uniform);
    }
    double R = 0;
    double PSI = 0;
    double time_elapsed = 0;
    for (size_t i = 0; i < iters; i++) {
        double dt = 1.0 / local_rates.sum;
        R += get_op(local_states) * dt;
        PSI += get_psi_op(local_states, local_rates) * dt;
        time_elapsed += dt;
        transition_site(local_states, local_deltas, local_rates, RNG, uniform);
    }
    Trial trial;
    trial.r = R / time_elapsed;
    trial.psi = PSI / time_elapsed;

    return trial;
}

Trial run_trial(
        size_t iters,
        size_t burn,
        States &local_states,
        Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG,
        Uniform &uniform
    ) {

    for (size_t i = 0; i < burn; i++) {
        transition_site(local_states, local_deltas, local_rates, RNG, uniform);
    }
    /* measure frequency of oscillations by counting how many times
       population 0 crosses the threshold. After detecting such a
       crossing we need to wait for a cooldown period to expire due to
       fluctuations near the crossing. */
    double period_start = 0;
    double period_end = 0;
    bool started = false;
    bool is_on_cooldown = false;
    int cd_timer = 0;
    int crossings = 0;
    size_t nprev = local_states.pop[0];
    const int cooldown = N * 1.8;
    const float threshold = N / 3.0;

    double R = 0;
    double PSI = 0;
    double omega;
    double time_elapsed = 0;
    for (size_t i = 0; i < iters; i++) {
        transition_site(local_states, local_deltas, local_rates, RNG, uniform);
        double dt = 1.0 / local_rates.sum;
        R += get_op(local_states) * dt;
        PSI += get_psi_op(local_states, local_rates) * dt;
        time_elapsed += dt;

        size_t n = local_states.pop[0];
        if (is_crossing(nprev, n, threshold, is_on_cooldown)) {
            if (!started) {
                started = true;
                period_start = time_elapsed;
            }
            period_end = time_elapsed;
            is_on_cooldown = true;
            crossings++;
        }
        if (is_on_cooldown) {
            cd_timer++;
        }
        if (cd_timer > cooldown) {
            is_on_cooldown = false;
            cd_timer = 0;
        }
        nprev = n;
    }
    if (crossings > 1) {
        omega = (crossings - 1) / 2.0 / (period_end - period_start);
    } else {
        omega = 0; // it is possible that no crossings ever happened
    }

    Trial trial;
    trial.r = R / time_elapsed;
    trial.psi = PSI / time_elapsed;
    trial.omega = omega;
    trial.duration = time_elapsed;

    return trial;
}

bool is_crossing(size_t nprev, size_t n, float t, bool is_on_cooldown) {
    if (((nprev <= t && n > t) || (nprev > t && n <= t)) && !is_on_cooldown) {
        return true;
    } else {
        return false;
    }
}

Batch run_batch(
            double coupling,
            size_t trial_iters,
            size_t trial_burn,
            size_t trials,
            bool verbose = false
        ) {
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &current_time);

    // shared variables for each trial
    const size_t seed = 23u * current_time.tv_nsec;
    Uniform uniform(0.0, 1.0);

    // private variables for each trial
    double r = 0;
    double r2 = 0;
    double psi = 0;
    double psi2 = 0;
    double omega = 0;
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    initialize_rates_table(coupling);
#pragma omp parallel default(none) \
    firstprivate(trial_iters,trial_burn,trials,uniform) \
    reduction(+:r,r2,psi,psi2,omega)
    {
        pcg32 RNG(seed);

        States states;
        Deltas deltas;
        Rates rates;
        reset_system(states, deltas, rates, RNG);
    #pragma omp for
        for (size_t i = 0; i < trials; i++) {
            pcg32 trial_rng(seed, i);
            initialize_states(states, RNG);
            initialize_deltas(states, deltas);
            initialize_rates(deltas, rates);
            Trial trial = run_trial(
                    trial_iters,
                    trial_burn,
                    states,
                    deltas,
                    rates,
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

std::string getDefaultTrialFilename(double coupling) {
    char fname[100];
    sprintf(
            fname,
            "trial-N-%05dK-%04dp-%6.6fa-%6.6f_v0.dat",
            N, K, p, coupling
        );
    fname[22] = '_';
    fname[32] = '_';
    int counter = 1;
    while (std::ifstream(fname)) {
        sprintf(
                fname,
                "trial-N-%05dK-%04dp-%6.6fa-%6.6f_v%d.dat",
                N, K, p, coupling, counter
            );
        fname[22] = '_';
        fname[32] = '_';
        counter++;
    }
    return fname;
}

std::string getDefaultBatchFilename(
        double coupling_start,
        double coupling_end
    ) {
    char fname[100];
    sprintf(
            fname,
            "batches-N-%05dK-%04dp-%6.6fa-%3.3f-%3.3f_v0.dat",
            N, K, p, coupling_start, coupling_end
        );
    fname[24] = '_';
    fname[34]  ='_';
    fname[40] = '_';
    int counter = 1;
    while (std::ifstream(fname)) {
        sprintf(
                fname,
                "batches-N-%05dK-%04dp-%6.6fa-%3.3f-%3.3f_v%d.dat",
                N, K, p, coupling_start, coupling_end, counter
            );
        fname[24] = '_';
        fname[34] = '_';
        fname[40] = '_';
        counter++;
    }
    return fname;
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
            if (cmdOptionExists(argv, argv+argc, "-ts")) {
                opt = getCmdOption(argv, argv+argc, "-ts");
                t_params.seed = stoi(opt);
            }
            if (cmdOptionExists(argv, argv+argc, "-tr")) {
                opt = getCmdOption(argv, argv+argc, "-tr");
                t_params.stream = stoi(opt);
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
                        b_params.coupling_start, b_params.coupling_end
                    );
            }
        } catch(...) {
            std::cerr << "Batch parameters not understood!\n"
                         "See -h for help.\n";
            return 0;
        }
    }

    std::cout << "Welcome, to Jurassick Park!\n";
    std::cout << "N=" << N << " K=" << K << " p=" << p
              << " topology_seed=" << TOPOLOGY_SEED << "\n\n";

    // RUN A BENCHMARK
    if (cmdOptionExists(argv, argv+argc, "--benchmark")) {
        struct timespec start, finish;
        double elapsed;
        double lowest = std::numeric_limits<double>::infinity();
        double highest = std::numeric_limits<double>::infinity();
        double avg = 0;
        int bench_trials = std::min(100, (int) (800000/(N*K) + 1));

        int iters = 10*N*log(N);
        States states;
        Deltas deltas;
        Rates rates;
        pcg32 RNG(DEFAULT_SEED< DEFAULT_STREAM);
        Uniform uniform(0.0, 1.0);
        initialize_everything(2.0, states, deltas, rates, RNG);
        for (int i=0; i<bench_trials; i++) {
            clock_gettime(CLOCK_MONOTONIC, &start);
            run_trial(iters, 0, states, deltas, rates, RNG, uniform);
            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            if (elapsed > highest) highest = elapsed;
            if (elapsed < lowest) lowest = elapsed;
            avg += elapsed;
        }
        avg /= bench_trials;

        std::cout << std::fixed << std::setprecision(6)
                  << "bench_trials=" << bench_trials << "\n"
                  << "Min:" << lowest << "s"
                  << " Max:" << highest << "s"
                  << " Avg:" << avg << "s" << "\n";
    }

    // RUN A TRIAL AND LOG IT TO A FILE
    if (cmdOptionExists(argv, argv+argc, "-t")) {

        States local_states;
        Deltas local_deltas;
        Rates local_rates;
        pcg32 RNG(t_params.seed, t_params.stream);
        Uniform uniform(0.0, 1.0);
        initialize_everything(
                t_params.coupling,
                local_states,
                local_deltas,
                local_rates,
                RNG
            );
        struct timespec start, finish;    // measure code run-times
        double elapsed;
        std::cout << "\nLogging trial to file. params:\n"
                  << "filename: " << t_params.filename << '\n'
                  << "coupling: " << t_params.coupling << '\n'
                  << "iters:    " << t_params.iters << '\n'
                  << "burn:     " << t_params.burn << '\n'
                  << "seed:     " << t_params.seed << '\n'
                  << "stream:   " << t_params.stream << '\n';

        FILE* trial_log_file = std::fopen(t_params.filename.c_str(), "w");
        fprintf(
                trial_log_file,
                "Graph_parameters: N=%d K=%d p=%f seed=%d\n"
                "Dynamics_parameters: coupling=%f iters=%lu burn=%lu "
                "seed=%lu stream=%lu\n"
                "r,psi,pop0,pop1,time_elapsed,dt\n",
                N, K, p, TOPOLOGY_SEED,
                t_params.coupling, t_params.iters, t_params.burn,
                t_params.seed, t_params.stream
            );
        clock_gettime(CLOCK_MONOTONIC, &start);
        log_trial_to_file(
                t_params.iters,
                t_params.burn,
                local_states,
                local_deltas,
                local_rates,
                RNG,
                uniform,
                trial_log_file
            );
        clock_gettime(CLOCK_MONOTONIC, &finish);
        std::fclose(trial_log_file);

        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        std::cout << "Trial finished in: " << elapsed << "s \n";
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

        // variables for measuring code run-times
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
