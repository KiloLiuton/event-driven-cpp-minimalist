#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#define SIN_PHI1 0.8660254037844387
#define DEFAULT_SEED 42u
#define DEFAULT_STREAM 23u

#include <random>
#include "default_topology.hpp"
#include "pcg_random/pcg_random.hpp"

typedef struct {
    uint8_t array[N];
    uint16_t pop[3];
} States;
typedef struct {
    double array[N];
    double sum;
} Rates;
typedef int16_t Deltas[N];
typedef std::uniform_real_distribution<double> Uniform;

// parameters to run a trial
struct trial_params {
    double coupling;
    size_t seed;
    size_t stream;
    size_t iters;
    size_t burn;
};

// data returned by a trial function
typedef struct {
    double r;
    double psi;
    double omega;
    double duration;
} Trial;

// parameters to run a batch
struct batch_params {
    double coupling_start;
    double coupling_end;
    int n_batches;
    int batch_id;
    size_t trials;
    size_t iters;
    size_t burn;
};

// data returned by a trial function
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

/* generates a new random configuration and update all dependencies */
void initialize_everything(
        double coupling,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG,
        std::string initial_condition,
        bool verbose
    );
/* reset the current lattice without changing the coupling strength */
void reset_system(
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG,
        std::string initial_condition
    );
/* populates rates table with a given coupling value */
void initialize_rates_table( double coupling, double rates_table[]);
/* get the transition rate based on the current delta value of site i */
double get_rate_from_table(uint16_t i, int16_t d, double rates_table[]);
/* populate a given states vector with a random configuration */
void initialize_random_states(States &local_states, pcg32 &RNG);
void initialize_uniform_states(States &local_states);
void initialize_wave_states(States &local_states, int num_bins);
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

/* run a trial with no logging of the omega parameter */
Trial run_no_omega_trial(
        size_t iters, size_t burn,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates,
        pcg32 &RNG, Uniform &uniform
    );

/* run a trial and return the average order parameters and frequency omega */
Trial run_trial(
        size_t iters, size_t burn,
        States &local_states, Deltas &local_deltas,
        Rates &local_rates, double rates_table[],
        pcg32 &RNG, Uniform &uniform
    );

/* determine if there is a crossing of threshold for frequency detection */
bool is_crossing(size_t nprev, size_t n, float t, bool is_on_cooldown);

/* execute a batch of trials and record the average order parameter */
Batch run_batch(
        double coupling,
        size_t trial_iters, size_t trial_burn, size_t trials,
        std::string initial_condition,
        bool verbose
    );

#endif
