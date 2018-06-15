#!/usr/bin/env python
from time import time
import networkx as nx
import random
import sys
import os


def writeHeader(fname, N, K, p, seed, INDEXES, NUMBER_OF_NEIGHBORS, NEIGHBOR_LIST):
    with open(fname, 'w') as f:
        K_MAX = max(NUMBER_OF_NEIGHBORS)
        K_MIN = min(NUMBER_OF_NEIGHBORS)
        f.write('#ifndef TOPOLOGY_H\n#define TOPOLOGY_H\n\n#include <iostream>\n\n')
        f.write('constexpr uint16_t N = %d\n;' % N)
        f.write('constexpr uint16_t K = %d\n;' % K)
        f.write('constexpr float p = %d\n;' % p)
        f.write('constexpr uint16_t K_MAX = %d;\n' % K_MAX)
        f.write('constexpr uint16_t K_MIN = %d;\n' % K_MIN)
        f.write('constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = %d;\n' % ((K_MAX - K_MIN + 1) * (K_MAX + K_MIN + 1)))
        f.write('constexpr uint32_t TOPOLOGY_SEED = %d;\n\n' % seed)
        f.write('constexpr uint32_t INDEXES[] = {\n')
        f.write(str(INDEXES).strip('[]') + '\n')
        f.write('};\n')
        f.write('constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {\n')
        f.write(str(NUMBER_OF_NEIGHBORS).strip('[]') + '\n')
        f.write('};\n')
        f.write('constexpr uint16_t NEIGHBOR_LIST[] = {\n')
        f.write(str(NEIGHBOR_LIST).strip('[]') + '\n')
        f.write('};\n')
        f.write('#endif')


def createHeader(N, K, p, s):
    random.seed(s)
    header_filename = ('%05d'%N + '-'
                     + '%04d'%K + '-'
                     + '_'.join(str(p).split('.')) + '-seed_'
                     + '%d'%s + '.h')

    if header_filename in os.listdir():
        print('File already exists! No need to create it :)')
        sys.exit()

    st = time()
    print('Ring Lattice...', end = ' ')
    G = nx.Graph()
    for i in range(N):
        for j in range(i - K, i + K + 1):
            if i == j: continue
            G.add_edge(i, j % N)
    en = time()
    print('Done! Took %.2f s' % (en - st))

    st = time()
    print('Rewiring...', end = ' ')
    if p > 0:
        S = set(range(N))
        for x in range(K):
            for i in range(N):
                if random.random() < p:
                    j = (i + x + 1) % N
                    k = random.choice(list(S - set(G[i]).union({i})))
                    G.remove_edge(i, j)
                    G.add_edge(i, k)
    en = time()
    print('Done! Took %.2f s' % (en - st))

    st = time()
    print('Parsing...', end=' ')
    NEIGHBOR_LIST = []
    INDEXES = [0]
    NUMBER_OF_NEIGHBORS = []
    cumsum = 0
    for i in range(N):
        nbrs = list(G[i])
        k = len(nbrs)
        cumsum += k
        NEIGHBOR_LIST += nbrs
        INDEXES += [cumsum]
        NUMBER_OF_NEIGHBORS += [k]
    en = time()
    print('Done! Took %.2f s\n' % (en - st))

    print('Writing to file %s' % header_filename)
    writeHeader(header_filename, N, K, p, s, INDEXES, NUMBER_OF_NEIGHBORS, NEIGHBOR_LIST)
    print('All done.\n')

if __name__ == "__main__":
    try:
        N, K, p = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
        s = 42
        if len(sys.argv) == 5:
            s = int(sys.argv[4])
        else:
            s = 42
    except:
        sys.exit('Usage - ./gen_header.py N K p s=42 (optional)')
    createHeader(N, K, p, s)
