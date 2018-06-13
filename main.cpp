#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <random>
#include <math.h>
#include <pcg_random.hpp>

// include the desired topology file and compile before running
#include "2000-200-0_0-seed_42.h"
// GLOBAL VARIABLES declared in topology header:
// uint16_t NEIGHBOR_LIST[];
// uint32_t INDEXES[];
// uint16_t NUMBER_OF_NEIGHBORS[];
// uint16_t K_MAX;
// uint16_t K_MIN;
// uint16_t N;
// uint16_t K;
// float p;
// uint32_t NUM_POSSIBLE_TRANSITIONS;

// GLOBAL VARIABLES
struct states {
    uint8_t array[N];
    uint16_t pop[3];
} states;
struct rates {
    double array[N];
    double sum;
} rates;
int16_t deltas[N];
double ratesTable[NUM_POSSIBLE_TRANSITIONS];
double timeElapsed;
std::uniform_real_distribution<double> UNIFORM(0,1);
pcg32 RNG(42u, 52u);

// INITIALIZER & GETTER FUNCTIONS
/* generates a new random configuration and update all dependencies */
void initialize_everything(double coupling);
/* populate the states vector with a random configuration */
void initialize_states();
/* populates delta vector (states must be populated) */
void initialize_deltas();
/* populates rates table */
void initialize_rates_table(double coupling);
/* get the transition rate for a site by accessing deltas vector and then
 * building the access index of the ratesTable. */
double get_rate_from_table(int site);
/* populates the rates vector */
void initialize_rates();
/* calculates de order parameter for the current state of the system */
double get_order_parameter();

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
/* execute a trial:
 * ITERS (pos integer) - execute ITERS steps with logging
 * BURN  (pos integer) - execute BURN steps without logging
 * reset (boolean)     - if false the current state of the lattice will be used
 *                       (this allows you to resume a simulation)
 * */
void run_trial(
    size_t ITERS,
    size_t BURN = 0,
    bool reset = true
);

// PRINTER FUNCTIONS
/* print the state for each site to stdout */
void print_states();
/* print deltas for all sites */
void print_deltas();
/* print transitions rates of every site */
void print_rates();
/* print a header to the trial output file (outputs results of a single run) */
void print_file_header(FILE* file, double coupling);

void initialize_everything(double a) {
    timeElapsed = 0;
    std::cout << "\nInitializing: "
              << "N=" << N
              << "  K=" << K
              << "  p=" << p
              << "  a=" << a << '\n';
    initialize_rates_table(a);
    std::cout << "Rates table initialized!\n";
    initialize_states();
    std::cout << "States initialized!\n";
    initialize_deltas();
    std::cout << "Deltas initialized!\n";
    initialize_rates();
    std::cout << "Rates initialized!\n\n";
}

void initialize_states() {
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
        const uint16_t ind = INDEXES[i];
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

double get_order_parameter() {
    uint16_t N0 = states.pop[0];
    uint16_t N1 = states.pop[1];
    uint16_t N2 = states.pop[2];
    return sqrt(N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / N;
}

void update_site(uint16_t i) {
    std::cout << "UPDATE SITE CALLED\n";
    uint8_t state = states.array[i];
    uint8_t nextState = (state + 1) % 3;
    uint16_t ki = NUMBER_OF_NEIGHBORS[i];
    states.array[i] = nextState;
    states.pop[state]--;
    states.pop[nextState]++;
    int16_t newDelta = 0;
    for (uint16_t j = INDEXES[i]; j < INDEXES[i] + ki; j++) {
        uint8_t nbState = states.array[NEIGHBOR_LIST[j]];
        if (nbState == nextState) {
            newDelta--;
        } else if (nbState == (nextState + 1) % 3) {
            newDelta++;
        }
    }
    deltas[i] = newDelta;
    rates.sum -= rates.array[i];
    rates.array[i] = ratesTable[(ki - K_MIN) * (ki + K_MAX) + ki + newDelta];
    rates.sum += rates.array[i];
}

void update_neighbors(uint16_t i) {
    std::cout << "UPDATE NEIGHBORS CALLED\n";
    // here we assume that site 'i' has already undergone transition
    uint8_t newState = states.array[i];
    for (uint16_t j = INDEXES[i]; j < INDEXES[i] + NUMBER_OF_NEIGHBORS[i]; j++) {
        uint16_t nbIndex = NEIGHBOR_LIST[j];
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
    std::cout << "TRANSITION INDEX CALLED\n";
    double partialRate = 0;
    double g = 0;
    double randomRate = UNIFORM(RNG) * rates.sum;
    for (uint16_t id = 0; id < N; ++id) {
        g = rates.array[id];
        partialRate += g;
        if (partialRate > randomRate) return id;
    }
    rates.sum = 0;
    for (uint16_t i = 0; i < N; i++) {
        rates.sum += rates.array[i];
    }
    std::cout << "Rates refreshed due to Overflow.";
    return N - 1;
}

void transition_site() {
    std::cout << "TRANSITION CALLED\n";
    uint16_t i = transitionIndex();
    timeElapsed += rates.sum;
    update_site(i);
    update_neighbors(i);
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

void print_file_header(FILE* f, double a) {
    fprintf(f, "N=%d k=%d p=%.6f a=%.6f\n", N, K, p, a);
}

void greeting() {
    std::cout << "Initial populations: ";
    std::cout << states.pop[0] << ' '
              << states.pop[1] << ' '
              << states.pop[2] << '\n';
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Initial order parameter: " << get_order_parameter() << '\n';
}

int main(int argc, char* argv[]) {
    // Simulation parameters
    const unsigned int random_stream = 2u;
    const unsigned int dynamics_seed = 20u;
    //const size_t BURN = 2*N*std::log(N);
    //const size_t ITERS = 2*N*std::log(N);
    const size_t BURN = 1;
    const size_t ITERS = 1;
    const size_t SAVE_INTERVAL = 1;
    const double coupling = 1.8;

    RNG.seed(dynamics_seed, random_stream);

    initialize_everything(coupling);

    greeting();
    std::cout << "Starting simulation with:"
              << "  BURN=" << BURN
              << "  ITERS=" << ITERS << '\n';

    for (size_t i = 0; i < BURN; ++i) {
        std::cout << "BURNING\n";
        transition_site();
    }

    char file_name[50];
    sprintf(file_name, "N-%05dK-%04dp-%3.3fa-%3.3f.dat", N, K, p, coupling);
    file_name[16] = '_';
    file_name[23] = '_';
    std::cout << "FOOOO" << file_name << '\n';

    FILE* OParameterLog;
    int counter = 1;
    while (std::ifstream(file_name)) {
        // file exists
        sprintf(file_name, "N-%05dK-%04dp-%3.3fa-%3.3f(%d).dat", N, K, p, coupling, counter);
        file_name[16] = '_';
        file_name[23] = '_';
        counter++;
    }
    OParameterLog = std::fopen(file_name, "w");
    print_file_header(OParameterLog, coupling);

    uint16_t log_counter = 1;
    for (size_t i = 0; i < ITERS; ++i) {
        transition_site();
        log_counter++;
        if (log_counter == SAVE_INTERVAL) {
            double r = get_order_parameter();
            std::fprintf(OParameterLog, "%6.6f,%d,%d\n", r, states.pop[0], states.pop[1]);
            log_counter = 1;
        }
    }

    std::fclose(OParameterLog);

    return 0;
}
