#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <algorithm>


uint8_t median(size_t n0, size_t n1, size_t n2)
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

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::string getCmdOption(char** begin, char** end, const std::string& option) {
    char** s = std::find(begin, end, option);
    if (s != end && ++s != end) {
        return *s;
    }
    return std::string();
}

std::string replaceAll(std::string& s, const std::string& a, const std::string& b) {
    size_t iter = s.find(a);
    while (iter != std::string::npos) {
        s.replace(iter, a.length(), b);
        iter = s.find(a);
    }
    return s;
}

#endif
