#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

#include <iostream>
#include <random>
#include <omp.h>
#include "pcg_random/pcg_random.hpp"
#include "dynamics.hpp"

// FUNCTION DEFINITIONS
void initialize_everything(
        double a,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG,
        std::string initial_condition="random",
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
    if (initial_condition == "random") {
        initialize_random_states(local_states, RNG);
    } else if (initial_condition == "uniform") {
        initialize_uniform_states(local_states, RNG);
    } else {
        initialize_random_states(local_states, RNG);
    }
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
        pcg32 &RNG,
        std::string initial_condition="random"
    ) {
    if (initial_condition == "random") {
        initialize_random_states(local_states, RNG);
    } else if (initial_condition == "uniform") {
        initialize_uniform_states(local_states, RNG);
    } else {
        initialize_random_states(local_states, RNG);
    }
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

void initialize_random_states(
        States &local_states,
        pcg32 &RNG
    ) {
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

void initialize_uniform_states(
        States &local_states,
        pcg32 &RNG
    ) {
        local_states.pop[0] = N;
        local_states.pop[1] = 0;
        local_states.pop[2] = 0;
        for (uint16_t i = 0; i < N; i++) {
            local_states.array[i] = 0;
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
        // TODO benchmark to see which is faster
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
    // TODO benchmark to see which is faster
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

Trial run_no_omega_trial(
        size_t iters, size_t burn,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG, Uniform &uniform,
        std::string initial_condition="random"
    ) {
    if (initial_condition == "uniform") {
        initialize_uniform_states(local_states, RNG);
    } else if (initial_condition == "random") {
        initialize_uniform_states(local_states, RNG);
    } else {
        initialize_random_states(local_states, RNG);
    }
    initialize_deltas(local_states, local_deltas);
    initialize_rates(local_deltas, local_rates, rates_table);
    for (size_t i = 0; i < burn; i++) {
        transition_site(
                local_states, local_deltas,
                local_rates, rates_table,
                RNG, uniform
            );
    }
    double R = 0;
    double PSI = 0;
    double time_elapsed = 0;
    for (size_t i = 0; i < iters; i++) {
        double dt = 1.0 / local_rates.sum;
        R += get_op(local_states) * dt;
        PSI += get_psi_op(local_states, local_rates) * dt;
        time_elapsed += dt;
        transition_site(
                local_states, local_deltas,
                local_rates, rates_table,
                RNG, uniform
            );
    }
    Trial trial;
    trial.r = R / time_elapsed;
    trial.psi = PSI / time_elapsed;

    return trial;
}

Trial run_trial(
        size_t iters, size_t burn,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG, Uniform &uniform
    ) {

    for (size_t i = 0; i < burn; i++) {
        transition_site(
                local_states, local_deltas,
                local_rates, rates_table,
                RNG, uniform
            );
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
        transition_site(
                local_states, local_deltas,
                local_rates, rates_table,
                RNG, uniform
            );
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
            size_t trial_iters, size_t trial_burn, size_t trials,
            std::string initial_condition="random",
            bool verbose = false
        ) {
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &current_time);

    // shared variables for each trial
    const size_t seed = 23u * current_time.tv_nsec;
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    Uniform uniform(0.0, 1.0);

    // private variables for each trial
    double r = 0;
    double r2 = 0;
    double psi = 0;
    double psi2 = 0;
    double omega = 0;
    initialize_rates_table(coupling, rates_table);
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    omp_set_num_threads(8);
#pragma omp parallel default(none) \
    shared(rates_table,initial_condition,std::cout) \
    firstprivate(trial_iters,trial_burn,trials,uniform,verbose) \
    reduction(+:r,r2,psi,psi2,omega)
    {
#pragma omp single
        {
            if (verbose) {
                int thrdnum = omp_get_num_threads();
                std::cout << "Batch running on " << thrdnum << " threads\n";
            }
        }
        pcg32 RNG(seed);

        States states;
        Deltas deltas;
        Rates rates;
        reset_system(states, deltas, rates, rates_table, RNG);
#pragma omp for
        for (size_t i = 0; i < trials; i++) {
            pcg32 trial_rng(seed, i);
            if (initial_condition == "random") {
                initialize_random_states(states, RNG);
            } else if (initial_condition == "uniform") {
                initialize_uniform_states(states, RNG);
            } else {
                initialize_random_states(states, RNG);
            }
            initialize_deltas(states, deltas);
            initialize_rates(deltas, rates, rates_table);
            Trial trial = run_trial(
                    trial_iters, trial_burn,
                    states, deltas,
                    rates, rates_table,
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
