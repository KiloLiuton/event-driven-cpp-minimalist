#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <algorithm>

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
