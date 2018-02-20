#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <math.h>
#include <pcg_random.hpp>

// 32 bits for indexing can go up to 2.5 billion connections
// topology constants: Regular Ring, eight neighbors per site
constexpr uint16_t NEIGHBOR_LIST[] = {
12,13,14,15,1,2,3,4,13,14,15,0,2,3,4,5,14,15,0,1,3,4,5,6,15,0,1,2,4,5,6,7,0,1,2,3,
5,6,7,8,1,2,3,4,6,7,8,9,2,3,4,5,7,8,9,10,3,4,5,6,8,9,10,11,4,5,6,7,9,10,11,12,5,6,
7,8,10,11,12,13,6,7,8,9,11,12,13,14,7,8,9,10,12,13,14,15,8,9,10,11,13,14,15,0,9,10,
11,12,14,15,0,1,10,11,12,13,15,0,1,2,11,12,13,14,0,1,2,3
};
constexpr uint32_t INDEXES[] = {
0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120
};
constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {
8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8
};

constexpr uint16_t K_MAX = 8; // largest element of NUMBER_OF_NEIGHBORS
constexpr uint16_t K_MIN = 8; // smallest element of NUMBER_OF_NEIGHBORS
constexpr uint16_t N = sizeof(INDEXES)/sizeof(INDEXES[0]);
constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = (K_MAX - K_MIN + 1)*
                                              (K_MAX + K_MIN + 1);

// lattice variables
struct states {
    uint8_t array[N];
    uint16_t N0, N1, N2;
} states;

int16_t deltas[N];

struct rates {
    double array[N];
    double sum;
} rates;

double ratesTable[NUM_POSSIBLE_TRANSITIONS];

// generates a new random congiguration and update all dependencies
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
// print the state for each site to stdout
void print_states();

void initialize_everything(double coupling_strength) {
    initialize_rates_table(coupling_strength);
    initialize_states();
    initialize_deltas();
    initialize_rates();
}

void initialize_states() {
    for (int i = 0; i < N; i++) {
        uint8_t state = rand() % 3;
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
        const int K = NUMBER_OF_NEIGHBORS[i];
        const int ind = INDEXES[i];
        const int8_t state = states.array[i];
        const int8_t nextState = (state+1)%3;
        int D = 0;
        for (int i = 0; i < K; i++) {
            int8_t nstate = states.array[NEIGHBOR_LIST[ind+i]];
            if (state == nstate) {
                D--;
            } else if (nextState == nstate) {
                D++;
            }
        }
        deltas[i] = D;
    }
}

void initialize_rates_table(double a) {
    int K = K_MIN;
    int D = -K_MIN;
    for (uint32_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        if (D > K) {
            K++;
            D = -K;
        }
        ratesTable[i] = std::exp(a*D/K);
        D++;
    }
}

double get_rate_from_table(int i) {
    const int K = NUMBER_OF_NEIGHBORS[i];
    const int D = deltas[i];
    int idx = (K - K_MIN) * (K + K_MIN) + K + D;
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
    std::cout << "System size: " << N << '\n';
    std::cout << "Initial states: ";
    print_states();
    std::cout << "Initial rates: ";
    print_rates();

    std::cout << "Initial OP: " << get_order_parameter() << '\n';

    pcg32 rng(42u, 52u);
    std::cout << "Dice rolls: " << '\n';
    for (int i = 0; i < 10; i++) {
        std::cout << rng(6) << ' ';
    }
    std::cout << '\n';

    return 0;
}
