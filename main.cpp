#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <random>
#include <math.h>
#include <pcg_random.hpp>

// include the desired topology file and compile before running
#include "7-2-0_0-seed_42.h"

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

std::uniform_real_distribution<double> uniform(0,1);

pcg32 rng(42u, 52u);

// INITIALIZER & GETTER FUNCTIONS
/* generates a new random configuration and update all dependencies */
void initialize_everything(double coupling_strength);
/* populate the states vector with a random configuration */
void initialize_states();
/* populates delta vector (states must be populated) */
void initialize_deltas();
/* populates rates table */
void initialize_rates_table(double coupling_strength);
/* get the transition rate for a site by accessing deltas vector and then
 * building the access index of the ratesTable. */
double get_rate_from_table(int site);
/* populates the rates vector */
void initialize_rates();
/* calculates de order parameter for the current state of the system */
double get_order_parameter();

// DYNAMICS FUNCTIONS
/* updates deltas, rates and states for site n */
void transition_site(int site_index);
/* update rates and delta values for all neighbors of a site which has been
 * updated */
void update_neighbors(int site_index);
/* select an index that will undergo transition */
int transitionIndex();

// PRINTER & DEBUGGING FUNCTIONS
/* print the state for each site to stdout */
void print_states();
/* print transitions rates of every site */
void print_rates();
/* print deltas for all sites */
void print_deltas();

void initialize_everything(double coupling_strength) {
    initialize_rates_table(coupling_strength);
    initialize_states();
    initialize_deltas();
    initialize_rates();
}

void initialize_states() {
    for (int i = 0; i < N; i++) {
        uint8_t state = rng(3);
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
    for (int i = 0; i < N; i++) {
        const int ki = NUMBER_OF_NEIGHBORS[i];
        const int ind = INDEXES[i];
        const int8_t state = states.array[i];
        const int8_t nextState = (state+1)%3;
        int d = 0;
        for (int j = 0; j < ki; j++) {
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
    int k = K_MIN;
    int d = -K_MIN;
    for (uint32_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        if (d > k) {
            k++;
            d = -k;
        }
        ratesTable[i] = std::exp(a*d/k);
        d++;
    }
}

double get_rate_from_table(int i) {
    const int k = NUMBER_OF_NEIGHBORS[i];
    const int d = deltas[i];
    int idx = (k - K_MIN) * (k + K_MIN) + k + d;
    return ratesTable[idx];
}

void initialize_rates() {
    rates.sum = 0;
    for (int i = 0; i < N; i++) {
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

void transition_site(int i) {
    int state = states.array[i];
    int nextState = (state + 1) % 3;
    uint16_t ki = NUMBER_OF_NEIGHBORS[i];
    states.array[i] = nextState;
    states.pop[state]--;
    states.pop[nextState]++;
    int newDelta = 0;
    for (size_t j = INDEXES[i]; j < INDEXES[i] + ki; j++) {
        int nbState = states.array[NEIGHBOR_LIST[j]];
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

void update_neighbors(int i) {
    // here we assume that site 'i' has already undergone transition
    int newState = states.array[i];
    for (size_t j = INDEXES[i]; j < INDEXES[i] + NUMBER_OF_NEIGHBORS[i]; j++) {
        int nbIndex = NEIGHBOR_LIST[j];
        int nbState = states.array[nbIndex];
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

int transitionIndex() {
    double partialRate = 0;
    double g = 0;
    double randomRate = uniform(rng) * rates.sum;
    for (int id = 0; id < N; ++id) {
        g = rates.array[id];
        partialRate += g;
        if (partialRate > randomRate) return id;
    }
    rates.sum = 0;
    for (int i = 0; i < N; i++) {
        rates.sum += rates.array[i];
    }
    return N - 1;
}

void print_states() {
    for (size_t i = 0; i < N; i++) {
        std::cout << unsigned(states.array[i]) << ' ';
    }
    std::cout << '\n';
}

void print_deltas() {
    for (size_t i = 0; i < N; i++) {
        std::cout << deltas[i] << ' ';
    }
    std::cout << '\n';
}

void print_rates() {
    std::cout << std::setprecision(2);
    for (size_t i = 0; i < N; i++) {
        std::cout << rates.array[i] << ' ';
    }
    std::cout << '\n';
}

int main(int argc, char* argv[]) {
    double coupling = 0.8;

    initialize_everything(coupling);

    // greeting message
    std::cout << "N = " << N << '\n';
    std::cout << "K = " << K << '\n';
    std::cout << std::fixed;
    std::cout << "Initial states: ";
    print_states();
    std::cout << "Initial deltas: ";
    print_deltas();
    std::cout << "Initial rates: ";
    print_rates();
    std::cout << "Initial OP: " << std::setprecision(4) << get_order_parameter() << '\n';

    transition_site(4);
    update_neighbors(4);
    std::cout << "\nUpdated states: ";
    print_states();
    std::cout << "Updated deltas: ";
    print_deltas();
    std::cout << "Updated rates: ";
    print_rates();
    std::cout << "Updated OP: " << std::setprecision(4) << get_order_parameter() << '\n';

    transition_site(4);
    update_neighbors(4);
    std::cout << "\nUpdated states: ";
    print_states();
    std::cout << "Updated deltas: ";
    print_deltas();
    std::cout << "Updated rates: ";
    print_rates();
    std::cout << "Updated OP: " << std::setprecision(4) << get_order_parameter() << '\n';

    transition_site(4);
    update_neighbors(4);
    std::cout << "\nUpdated states: ";
    print_states();
    std::cout << "Updated deltas: ";
    print_deltas();
    std::cout << "Updated rates: ";
    print_rates();
    std::cout << "Updated OP: " << std::setprecision(4) << get_order_parameter() << '\n';

    return 0;
}
