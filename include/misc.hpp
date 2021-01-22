#ifndef MISC_HPP
#define MISC_HPP

#define PI 3.14159265359

#include <string>
#include <algorithm>
#include <dynamics.hpp>
#include <pcg_random.hpp>


inline uint8_t median(size_t n0, size_t n1, size_t n2)
{
    // return the phase value of the majority, selecting randomly on ties.
    if (n0 > n2 && n0 > n1)       return 0;
    if (n1 > n2 && n1 > n0)       return 1;
    if (n2 > n0 && n2 > n1)       return 2;
    if ((n0 == n1) && (n1 == n2)) return rand() % 3;
    if (n0 == n1)                 return rand() % 2;
    if (n1 == n2)                 return rand() % 2 + 1;
    if (n0 == n2)                 return (rand() % 2)? 0 : 2;
    return 0;
}

inline bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

inline std::string getCmdOption(char** begin, char** end, const std::string& option) {
    char** s = std::find(begin, end, option);
    if (s != end && ++s != end) {
        return *s;
    }
    return std::string();
}

inline std::string replaceAll(std::string& s, const std::string& a, const std::string& b) {
    size_t iter = s.find(a);
    while (iter != std::string::npos) {
        s.replace(iter, a.length(), b);
        iter = s.find(a);
    }
    return s;
}

inline void initialize_custom_natural_frequencies(NaturalFreqs &g, std::string nfreq_bias, double nfreq_bias_amplitude, pcg32 &RNG)
{
    if (nfreq_bias == "sinesqr") {
        for (int n=0; n<N; n++) {
            g[n] = 1.0 + 0.5 * (1.0 + sin(n*PI/(2*N))*sin(n*PI/(2*N)));
        }
    } else if (nfreq_bias == "linear") {
        for (int n=0; n<N/2; n++ ) { g[n] = 1.0 + 0.5 * (1.0 + (double) n/N); }
        for (int n=N/2; n>=0; n--) { g[n] = 1.0 + 0.5 * (1.0 - (double) n/N); }
    } else {
        initialize_natural_frequencies(1.0, nfreq_bias_amplitude, g, RNG);
    }
}

#endif
