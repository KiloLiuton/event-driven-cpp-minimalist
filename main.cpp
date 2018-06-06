#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <vector>
#include <random>
#include <math.h>
#include <pcg_random.hpp>
#include "topology.h"

// GLOBAL VARIABLES declared in topology.h
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
    uint16_t N0, N1, N2;
} states;

struct rates {
    double array[N];
    double sum;
} rates;

int16_t deltas[N];

double ratesTable[NUM_POSSIBLE_TRANSITIONS];

std::uniform_real_distribution<double> uniform(0,1);

pcg32 rng(42u, 52u);

// FUNCTIONS
// generates a new random configuration and update all dependencies
void initialize_everything(double coupling_strength);
// populate the states vector with a random configuration
void initialize_states();
// populates delta vector (states must be populated)
void initialize_deltas();
// populates rates table
void initialize_rates_table(double coupling_strength);
// get the transition rate for a site by accessing deltas vector and rateTable
double get_rate_from_table(int site);
// populates the rates vector
void initialize_rates();
// calculates de order parameter for the current state of the system
double get_order_parameter();
// select an index that will undergo transition
int transitionIndex();
// print the state for each site to stdout
void print_states();
// print transitions rates of every site
void print_rates();

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
            states.N0++;
            break;
        case 1:
            states.N1++;
            break;
        case 2:
            states.N2++;
            break;
        default:
            std::cout << "Error: Generated state out of range\n";
            break;
        }
    }
}

void initialize_deltas() {
    for (int i = 0; i < N; i++) {
        const int k = NUMBER_OF_NEIGHBORS[i];
        const int ind = INDEXES[i];
        const int8_t state = states.array[i];
        const int8_t nextState = (state+1)%3;
        int d = 0;
        for (int i = 0; i < k; i++) {
            int8_t nstate = states.array[NEIGHBOR_LIST[ind+i]];
            if (state == nstate) {
                d--;
            } else if (nextState == nstate) {
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
    uint16_t N0 = states.N0;
    uint16_t N1 = states.N1;
    uint16_t N2 = states.N2;
    return sqrt(N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / N;
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
    return N - 1;
}

void print_states() {
    for (size_t i = 0; i < N; i++) {
        std::cout << unsigned(states.array[i]) << ' ';
    }
    std::cout << '\n';
}

void print_rates() {
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
    std::cout << "Initial states: ";
    print_states();
    std::cout << "Initial rates: ";
    print_rates();
    std::cout << "Initial OP: " << get_order_parameter() << '\n';

    return 0;
}
