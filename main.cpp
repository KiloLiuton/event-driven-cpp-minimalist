#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <math.h>
#include <pcg_random.hpp>

// 32 bits for indexing can go up to ~ N = 50000 with n = 50000 (complete graph)

constexpr uint16_t neighbor_list[] = { 6,7,1,2,7,0,2,3,0,1,3,4,1,2,4,5,2,3,5,6,3,4,6,7,4,5,7,0,5,6,0,1 };
constexpr uint32_t indexes[] = { 0,4,8,12,16,20,24,28 };
constexpr uint16_t number_of_neighbors[] = { 4,4,4,4,4,4,4,4 };
constexpr uint16_t maxK = 4; // largest element of number_of_neighbors
constexpr uint16_t minK = 4; // smallest element of number_of_neighbors
constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = (2*minK+1)*(maxK-minK) +
                                              (maxK-minK)*(1+maxK-minK);

constexpr int N = sizeof(indexes)/sizeof(indexes[0]);
constexpr int n = number_of_neighbors[0];

uint8_t states[N];
uint16_t N0 = 0, N1 = 0, N2 = 0;
double rates[N];
double ratesTable[NUM_POSSIBLE_TRANSITIONS];

void initialize_states()
{
    for (int i = 0; i < N; i++) {
        uint8_t state = rand() % 3;
        states[i] = state;
        switch (state)
        {
        case 0:
            N0++;
            break;
        case 1:
            N1++;
            break;
        case 2:
            N2++;
            break;
        default:
            std::cout << "Error: Generated state out of range\n";
            break;
        }
    }
}

void initialize_rates_table(double a)
{
    int K = minK;
    int D = -minK;
    for (uint32_t i = 0; i < NUM_POSSIBLE_TRANSITIONS; i++) {
        if (D > K) {
            K++;
            D = -K;
        }
        ratesTable[i] = std::exp(a*D/K);
        D++;
    }
}

double compute_rate(int i)
{
    // compute the transition rate for site i
    double rate = 0;
    int D = 0;
    int K = number_of_neighbors[i];
    int ind = indexes[i];
    int8_t state = states[i];
    int8_t nextState = (state+1)%3;
    for (int i = 0; i < K; i++) {
        int8_t nstate = states[ neighbor_list[ind+i] ];
        if (state == nstate) {
            D--;
        } else if (nextState == nstate) {
            D++;
        }
    }
    // TODO: get rate for site i
    return rate;
}

void initialize_rates()
{
    for (int i = 0; i < N; i++) {
        rates[i] = compute_rate(i);
        continue;
    }
}

double get_order_parameter()
{
    return sqrt(N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2) / N;
}

void print_states()
{
    for (int i = 0; i < N; i++) {
        std::cout << unsigned(states[i]) << ' ';
    }
    std::cout << '\n';
}

int main(int argc, char* argv[])
{
    initialize_states();
    print_states();
    std::cout << get_order_parameter() << '\n';

    pcg32 rng(42u, 52u);
    std::cout << "Dice rolls: " << '\n';
    for (int i = 0; i < 10; i++) {
        std::cout << rng(6) << ' ';
    }
    std::cout << '\n';

    return 0;
}
