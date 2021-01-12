#include "dynamics.hpp"

class Benchmark {
private:
    double _fastest, _slowest, _average;
    int _iters;
public:
    Benchmark(int iters);
    void run();
    double fastest() { return _fastest; }
    double slowest() { return _slowest; }
    double average() { return _average; }
};

Benchmark::Benchmark(int iters) {
    _iters = iters;
    _fastest = std::numeric_limits<double>::infinity();
    _slowest = -std::numeric_limits<double>::infinity();
    _average = 0;
}

void Benchmark::run() {
    struct timespec start, finish;
    double elapsed;
    States states;
    Deltas deltas;
    Rates rates;
    double rates_table[NUM_POSSIBLE_TRANSITIONS];
    NaturalFreqs g;
    pcg32 RNG(DEFAULT_SEED, DEFAULT_STREAM);
    Uniform uniform(0.0, 1.0);
    initialize_everything(2.0, states, deltas, rates, rates_table, g, RNG, "random", false);
    for (int i=0; i<_iters; i++) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        run_trial(
                10*N*log(N), 0,
                states, deltas,
                rates, rates_table,
                g,
                RNG, uniform
            );
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        if (elapsed > _slowest)_slowest = elapsed;
        if (elapsed < _fastest) _fastest = elapsed;
        _average += elapsed;
    }
    _average /= _iters;
}
