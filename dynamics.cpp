#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

#include <pcg_random.hpp>
#include <random>

typedef struct {
    uint8_t array[N];
    uint16_t pop[3];
} States;
typedef int16_t Deltas[N];
typedef struct {
    double array[N];
    double sum;
} Rates;
typedef std::uniform_real_distribution<double> Uniform;

/* generates a new random configuration and update all dependencies */
void initialize_everything(
        double coupling,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG,
        bool verbose
    );
/* reset the current lattice without changing the coupling strength */
void reset_system(
        States &local_states, Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG
    );
/* populates rates table with a given coupling value */
void initialize_rates_table( double coupling, double rates_table[]);
/* get the transition rate based on the current delta value of site i */
double get_rate_from_table(uint16_t i, int16_t d, double rates_table[]);
/* populate a given states vector with a random configuration */
void initialize_states(States &local_states, pcg32 &RNG);
/* populates a given rates vector from the current rates table */
void initialize_rates(
        Deltas &local_deltas,
        Rates &local_rates, double rates_table[]
    );
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
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[]
    );
/* select an index that will undergo transition */
uint16_t transitionIndex(
        Rates &local_rates,
        pcg32 &RNG, Uniform &uniform
    );
/* performs a complete step of the event driven simulation */
uint16_t transition_site(
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG, Uniform &uniform
    );

// FUNCTION DEFINITIONS
void initialize_everything(
        double a,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
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
    initialize_rates_table(a, rates_table);
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
    initialize_rates(local_deltas, local_rates, rates_table);
    if (verbose) {
        std::cout << "Rates initialized!\n\n";
    }
}

void reset_system(
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG
    ) {
    initialize_states(local_states, RNG);
    initialize_deltas(local_states, local_deltas);
    initialize_rates(local_deltas, local_rates, rates_table);
}

void initialize_rates_table( double a, double rates_table[]) {
    uint16_t k = K_MIN;
    int16_t d = -K_MIN;
    for (uint32_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        if (d > k) {
            k++;
            d = -k;
        }
        rates_table[i] = std::exp(a*d/k);
        d++;
    }
}

double get_rate_from_table(uint16_t i, int16_t d, double rates_table[]) {
    const uint16_t k = NUMBER_OF_NEIGHBORS[i];
    uint32_t idx = (k - K_MIN) * (k + K_MIN) + k + d;
    return rates_table[idx];
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

void initialize_rates(
        Deltas &local_deltas,
        Rates &local_rates, double rates_table[]
    ) {
    local_rates.sum = .0;
    for (uint16_t i = 0; i < N; i++) {
        int d = local_deltas[i];
        const double r = get_rate_from_table(i, d, rates_table);
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
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[]
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
        double nbRate = get_rate_from_table(nbIndex, local_deltas[nbIndex], rates_table);
        rateIncrease += nbRate;
        rateDecrease -= local_rates.array[nbIndex];
        local_rates.array[nbIndex] = nbRate;
    }
    local_deltas[i] = siteDelta;
    double siteRate = get_rate_from_table(i, local_deltas[i], rates_table);
    rateIncrease += siteRate;
    rateDecrease -= local_rates.array[i];
    local_rates.array[i] = siteRate;
    local_rates.sum += rateIncrease + rateDecrease;
}

uint16_t transitionIndex(
        Rates &local_rates,
        pcg32 &RNG, Uniform &uniform
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
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG, Uniform &uniform
    ) {
    uint16_t i = transitionIndex(local_rates, RNG, uniform);
    update_site(i, local_states, local_deltas, local_rates, rates_table);
    return i;
}

#endif
