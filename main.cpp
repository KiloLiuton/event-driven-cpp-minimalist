#include <iostream>
#include <array>
#include <vector>
#include <math.h>

constexpr int N = 10;
constexpr int K = 2;

// 32 bits for indexing can go up to ~ N = 50000 with n = 50000 (complete graph)

uint8_t states[ N ];
uint16_t N0 = 0, N1 = 0, N2 = 0;
constexpr uint16_t neighbor_list[] = { 6,7,1,2,7,0,2,3,0,1,3,4,1,2,4,5,2,3,5,6,3,4,6,7,4,5,7,0,5,6,0,1 };
constexpr uint32_t indexes[] = { 0,4,8,12,16,20,24,28 };
constexpr uint16_t kernel_sizes[] = { 4,4,4,4,4,4,4,4 };

void initialize_states()
{
    for ( int i = 0; i < N; i++ )
    {
        uint8_t state = rand() % 3;
        states[ i ] = state;
        switch ( state )
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

double get_order_parameter()
{
    return sqrt( N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2 ) / N;
}

void print_states()
{
    for ( int i = 0; i < N; i++ )
    {
        std::cout << states[i] << ' ';
    }
    std::cout << '\n';
}

int main( int argc, char* argv[] )
{
    initialize_states();
    print_states();

    return 0;
}
