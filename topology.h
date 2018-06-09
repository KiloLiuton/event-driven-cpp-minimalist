#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <iostream>

constexpr uint16_t N = 16;
constexpr uint16_t K = 4; // neighbors to each side on the regular ring
constexpr float p = .0; // rewire probability

// 32 bits for indexing can go up to 2.5 billion connections
// topology constants: Regular Ring, eight neighbors per site
constexpr uint16_t NEIGHBOR_LIST[] = {
12,13,14,15,1,2,3,4,13,14,15,0,2,3,4,5,14,15,0,1,3,4,5,6,15,0,1,2,4,5,6,7,0,1,2,
3,5,6,7,8,1,2,3,4,6,7,8,9,2,3,4,5,7,8,9,10,3,4,5,6,8,9,10,11,4,5,6,7,9,10,11,12,
5,6,7,8,10,11,12,13,6,7,8,9,11,12,13,14,7,8,9,10,12,13,14,15,8,9,10,11,13,14,15,
0,9,10,11,12,14,15,0,1,10,11,12,13,15,0,1,2,11,12,13,14,0,1,2,3
};
constexpr uint32_t INDEXES[] = {
0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120
};
constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {
8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8
};

constexpr uint16_t K_MAX = 8; // largest element of NUMBER_OF_NEIGHBORS
constexpr uint16_t K_MIN = 8; // smallest element of NUMBER_OF_NEIGHBORS
constexpr uint32_t NUM_POSSIBLE_TRANSITIONS =
	(K_MAX - K_MIN + 1) * (K_MAX + K_MIN + 1);

#endif
