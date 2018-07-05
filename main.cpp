#define SIN_PHI1 0.8660254037844387
//#ifndef HEADER
//#define HEADER "20-3-0_0-seed_42.h"
//#endif

#include <iostream>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <random>
#include <math.h>
#include <pcg_random.hpp>
// #include "greeting.h"

// GLOBAL VARIABLES
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
    double omega = 0;
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
    size_t stream,
    bool verbose = false
);
/* populate the states vector with a random configuration */
void initialize_states();
/* populates delta vector (states must be populated) */
void initialize_deltas();
/* populates rates table (only necessary if coupling changed) */
void initialize_rates_table(double coupling);
/* get the transition rate for a site by accessing deltas vector and then
 * building the access index of the ratesTable. */
double get_rate_from_table(int site);
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
void update_neighbors(uint16_t site_index);
/* select an index that will undergo transition */
uint16_t transitionIndex();
/* performs a complete step of the event driven simulation */
void transition_site();
/* run a trial for ITERS time steps after burning BURN steps and save results
 * to log_file. Logged columns are:
 * r**2,N0,N1,time */
void log_trial_to_file(
        size_t ITERS,
        size_t BURN,
        size_t seed,
        size_t stream,
        FILE* log_file,
        size_t SAVE_INTERVAL
    );
/* run a trial and return the average order parameters after BURN iters */
trial run_trial(
        size_t ITERS,
        size_t BURN,
        size_t seed,
        size_t stream
    );
/* run a trial and return the average order parameters and frequency omega
 * after BURN iters */
trial run_omega_trial(
        size_t ITERS,
        size_t BURN,
        size_t seed,
        size_t stream
    );
/* execute a batch of trials and record the average order parameter */
batch run_batch(
        size_t TRIALS,
        size_t ITERS,
        size_t BURN,
        double a,
        bool RUN_OMEGA
    );

// PRINTER FUNCTIONS
/* print the state for each site to stdout */
void print_states();
/* print deltas for all sites */
void print_deltas();
/* print transitions rates of every site */
void print_rates();
/* print a header to the trial output file (outputs results of a single run) */
void print_file_header(FILE* file, double coupling, int BURN, int ITERS);
/* print a header to the batches output file (many trials for many couplings) */
void print_batches_file_header(FILE* file, int TRIALS, int BURN, int ITERS);

void initialize_everything(double a, size_t s1, size_t s2, bool verbose) {
    RNG.seed(s1, s2);
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
        int16_t d = 0;
        for (uint16_t j = 0; j < ki; j++) {
            int8_t nbState = states.array[NEIGHBOR_LIST[ind+j]];
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
    // TODO remove
    if (idx < 0 || idx >= NUM_POSSIBLE_TRANSITIONS) std::cout << "get_rate_from_table ALERT!!! BUG ALERT!!!\n";
    return ratesTable[idx];
}

void initialize_rates() {
    rates.sum = 0.0;
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
    uint16_t ki = NUMBER_OF_NEIGHBORS[i];
    states.array[i] = nextState;
    states.pop[state]--;
    states.pop[nextState]++;
    int16_t newDelta = 0;
    for (uint32_t j = INDEXES[i]; j < INDEXES[i] + ki; j++) {
        uint8_t nbState = states.array[NEIGHBOR_LIST[j]];
        if (nbState == nextState) {
            newDelta--;
        } else if (nbState == (nextState + 1) % 3) {
            newDelta++;
        }
    }
    // TODO remove
    std::cout << "atualizando sitio " << i << "\n";
    deltas[i] = newDelta;
    double oldrate = rates.array[i];
    // TODO remove
    if (std::isnan(oldrate)) std::cout << "oldrate foodeu\n";
    rates.sum -= oldrate;
    double newrate = ratesTable[(ki - K_MIN) * (ki + K_MAX) + ki + newDelta];

    const uint16_t k = NUMBER_OF_NEIGHBORS[i];
    const int16_t d = deltas[i];
    uint32_t idx = (k - K_MIN) * (k + K_MIN) + k + d;

    // TODO remove
    if (std::isnan(newrate)) std::cout << "newrate foodeu\n";
    rates.array[i] = newrate;
    rates.sum += rates.array[i];
    // TODO remove
    if (std::isnan(rates.sum)) std::cout << "Deu pau no sitio " << i << "\n";
}

void update_neighbors(uint16_t i) {
    // here we assume that site 'i' has already undergone transition
    uint8_t newState = states.array[i];
    for (uint32_t j = INDEXES[i]; j<INDEXES[i] + NUMBER_OF_NEIGHBORS[i]; j++) {
        uint16_t nbIndex = NEIGHBOR_LIST[j];
        // TODO remove
        if (nbIndex < 0 || nbIndex >= N) std::cout << "ALERT!!! BUG ALERT!!!\n";
        uint16_t nbState = states.array[nbIndex];
        if (nbState == newState) {
            deltas[nbIndex] -= 1;
        } else if (nbState == (newState + 1) % 3) {
            deltas[nbIndex] -= 1;
        } else {
            deltas[nbIndex] += 2;
        }
        // be sure to update rates only after updating `deltas` because
        // get_rate_from_table uses it to get the updated value for the rate
        rates.sum -= rates.array[nbIndex];
        rates.array[nbIndex] = get_rate_from_table(nbIndex);
        rates.sum += rates.array[nbIndex];
    }
}

uint16_t transitionIndex() {
    double partialRate = 0;
    double g = 0;
    double rn = UNIFORM(RNG);
    double randomRate = rn * rates.sum;
    // TODO remove
    std::cout << "FOO rn=" << rn << " randomRate=" << randomRate << " sum=" << rates.sum << '\n';
    for (uint16_t id = 0; id < N; ++id) {
        g = rates.array[id];
        partialRate += g;
        if (partialRate > randomRate) return id;
    }
    rates.sum = 0;
    for (uint16_t i = 0; i < N; i++) {
        rates.sum += rates.array[i];
    }
    std::cout << "Rates refreshed due to Overflow.\n";
    return N - 1;
}

void transition_site() {
    uint16_t i = transitionIndex();
    update_site(i);
    update_neighbors(i);
}

void log_trial_to_file(
            size_t I, size_t B,
            size_t seed, size_t stream,
            FILE* f, size_t SAVE_INTERVAL
        ) {
    std::cout << "Logging trial to file with:"
              << "  BURN=" << B
              << "  ITERS=" << I << '\n';

    size_t PROGRESS_INTERVAL = (I+B)/20;
    size_t log_counter = 1;
    size_t progress_counter = 1;
    double time_elapsed = 0;
    for (size_t i = 0; i < (I + B); ++i) {
        double dt = 1.0 / rates.sum;
        time_elapsed += dt;
        transition_site();

        if (log_counter == SAVE_INTERVAL) {
            double r = get_squared_op();
            std::fprintf(
                    f, "%16.16f,%d,%d,%f,%f\n", r, states.pop[0], states.pop[1], time_elapsed, dt
                );
            log_counter = 0;
        }
        if (progress_counter == PROGRESS_INTERVAL) {
            std::cout << std::fixed << std::setprecision(1)
                      << "[" << (float) i/(I+B)*100 << "%]\n";
            progress_counter = 0;
        }
        log_counter++;
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

trial run_omega_trial(
            size_t I, size_t B,
            size_t seed, size_t stream
        ) {
    RNG.seed(seed, stream);
    initialize_states();
    initialize_deltas();
    initialize_rates();
    for (size_t i = 0; i < B; i++) transition_site();

    const int cooldown = N * 1.8;
    float t_start = 0;
    float t_end = 0;
    bool started = false;
    bool is_on_cooldown = false;
    int crossings = 0;
    int counter = 0;
    size_t nprev = states.pop[0];
    const float thold = N / 3.0;

    double R = 0;
    double PSI = 0;
    double time_elapsed = 0;
    for (size_t i = 0; i < I; i++) {
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
    }
    return Trial;
}

batch run_batch(size_t TRIALS, size_t TRIAL_I, size_t TRIAL_B,
        double a, bool RUN_OMEGA) {
    initialize_everything(a, 42u, 52u); // seeds for the initial state
    size_t PROGRESS_INTERVAL = TRIALS / 5;
    size_t progress_counter = 1;
    double r = 0;
    double r2 = 0;
    double psi = 0;
    double omega = 0;
    for (size_t i = 0; i < TRIALS; i++) {
        size_t seed = i;
        size_t stream = i;
        trial Trial;
        if (RUN_OMEGA) {
            Trial = run_omega_trial(TRIAL_I, TRIAL_B, seed, stream);
        } else {
            Trial = run_trial(TRIAL_I, TRIAL_B, seed, stream);
        }
        r += Trial.r;
        r2 += std::pow(Trial.r, 2.0);
        psi += Trial.psi;
        omega += Trial.omega;

        if (progress_counter == PROGRESS_INTERVAL) {
            std::cout << std::setprecision(1) << std::fixed
                      << "[" << (float) i / TRIALS * 100 << "%]\n";
            progress_counter = 0;
        }
        progress_counter++;
    }
    batch Batch;
    Batch.r = r / TRIALS;
    Batch.r2 = r2 / TRIALS;
    Batch.psi = psi / TRIALS;
    Batch.omega = omega / TRIALS;
    return Batch;
}

void print_states() {
    for (uint16_t i = 0; i < N; i++) {
        std::cout << unsigned(states.array[i]) << ' ';
    }
    std::cout << '\n';
}

void print_deltas() {
    for (uint16_t i = 0; i < N; i++) {
        std::cout << deltas[i] << ' ';
    }
    std::cout << '\n';
}

void print_rates() {
    std::cout << std::setprecision(2);
    for (uint16_t i = 0; i < N; i++) {
        std::cout << rates.array[i] << ' ';
    }
    std::cout << '\n';
}

void print_file_header(FILE* f, double a, int BURN, int ITERS) {
    fprintf(f, "N=%d,k=%d,p=%.6f,a=%.6f,BURN=%d,ITERS=%d\n", N, K, p, a, BURN, ITERS);
    fprintf(f, "r**2,N0,N1,time,dt\n");
}

void print_batches_file_header(FILE* f, int TRIALS, int BURN, int ITERS) {
    fprintf(f, "N=%d, K=%d, p=%f, TRIALS=%d, BURN=%d, ITERS=%d\n", N, K, p, TRIALS, BURN, ITERS);
    fprintf(f, "<<r>>,<<r>^2>,<<psi>>,a\n");
}

void welcome(int argc, char* argv[],
        bool RUN_TRIAL, bool RUN_BATCH, bool RUN_OMEGA) {
    // std::cout << greeting;
    std::cout << "Welcome to the event driven simulation!\n";
    std::cout << "Parameters: ";
    std::cout << "N=" << N << "    K=" << K << "    p=" << p << '\n';
    std::cout << "K_MIN=" << K_MIN << "    K_MAX=" << K_MAX << '\n';
    std::cout << "Command line arguments: ";
    for (int i = 1; i < argc; i++) {
        std::cout << argv[i] << "    ";
    }
    std::cout << '\n';
    std::cout << "Selected modes:\n";
    std::cout << "RUN_TRIAL: " << (RUN_TRIAL ? "True":"False") << '\n';
    std::cout << "RUN_BATCH: " << (RUN_BATCH ? "True":"False") << '\n';
    std::cout << "RUN_OMEGA: " << (RUN_OMEGA ? "True":"False") << '\n';
}

bool is_arg(const char* arg, int argc, char* argv[]) {
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], arg) == 0) return true;
    }
    return false;
}

int main(int argc, char* argv[]) {
    welcome(argc, argv, false, false, false);
    const unsigned int seed = 20u;
    const unsigned int stream = 2u;
    initialize_everything(2.0, seed, stream, true);
    for (size_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        std::cout << "rt[" << i << "] = " << ratesTable[i] << '\n';
    }
    for (int i = 0; i < 10; i++) {
        transition_site();
    }

    return 0;
}
