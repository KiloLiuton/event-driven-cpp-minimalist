#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <iostream>

constexpr uint16_t N = 20;
constexpr uint16_t K = 3;
constexpr float p = 0.000000;
constexpr uint16_t K_MAX = 6;
constexpr uint16_t K_MIN = 6;
constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = 13;
constexpr uint32_t TOPOLOGY_SEED = 42;

constexpr uint32_t INDEXES[] = {
0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114
};
constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6
};
constexpr uint16_t NEIGHBOR_LIST[] = {
17, 18, 19, 1, 2, 3, 0, 18, 19, 2, 3, 4, 0, 1, 19, 3, 4, 5, 0, 1, 2, 4, 5, 6, 1, 2, 3, 5, 6, 7, 2, 3, 4, 6, 7, 8, 3, 4, 5, 7, 8, 9, 4, 5, 6, 8, 9, 10, 5, 6, 7, 9, 10, 11, 6, 7, 8, 10, 11, 12, 7, 8, 9, 11, 12, 13, 8, 9, 10, 12, 13, 14, 9, 10, 11, 13, 14, 15, 10, 11, 12, 14, 15, 16, 11, 12, 13, 15, 16, 17, 12, 13, 14, 16, 17, 18, 13, 14, 15, 17, 18, 19, 0, 14, 15, 16, 18, 19, 0, 1, 15, 16, 17, 19, 0, 1, 2, 16, 17, 18
};
#endif