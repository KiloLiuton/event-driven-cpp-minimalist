import networkx as nx
import random
import sys

random.seed(42)

try:
    N, K, p = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
except:
    sys.exit('Usage - ./gen_header.py N K p')

G = nx.Graph()

for i in range(N):
    for j in range(i - K, i + K + 1):
        if i == j: continue
        G.add_edge(i, j % N)

if p > 0:
    S = set(range(N))
    for x in range(K):
        for i in range(N):
            if random.random() < p:
                j = (i + x + 1) % N
                k = random.choice(list(S - set(G[i]).union({i})))
                G.remove_edge(i, j)
                G.add_edge(i, k)

NEIGHBOR_LIST = []
INDEXES = [0]
NUMBER_OF_NEIGHBORS = []
cumsum = 0
for n, nbrs in G.adj.items():
    nbrs = list(nbrs)
    k = len(nbrs)
    cumsum += k
    NEIGHBOR_LIST += nbrs
    INDEXES += [cumsum]
    NUMBER_OF_NEIGHBORS += [k]
print('DONE!')

def writeHeader():
    with open('header.h', 'w') as f:
        f.write('#ifndef TOPOLOGY_H\n#define TOPOLOGY_H\n\n#include <iostream>\n\n')
        f.write('constexpr uint16_t N = %d\n;' % N)
        f.write('constexpr uint16_t K = %d\n;' % K)
        f.write('constexpr float p = %d\n;' % p)
        f.write('constexpr uint16_t NEIGHBOR_LIST[] = {\n')
        f.write(str(NEIGHBOR_LIST).strip('[]') + '\n')
        f.write('};\n')
        f.write('constexpr uint16_t INDEXES[] = {\n')
        f.write(str(INDEXES).strip('[]') + '\n')
        f.write('};\n')
        f.write('constexpr uint16_t NUMBER_OF_NEIGHBORS[] = {\n')
        f.write(str(NUMBER_OF_NEIGHBORS).strip('[]') + '\n')
        f.write('};\n')
        K_MAX = max(NUMBER_OF_NEIGHBORS)
        K_MIN = min(NUMBER_OF_NEIGHBORS)
        f.write('constexpr uint16_t K_MAX = %d;\n' % K_MAX)
        f.write('constexpr uint16_t K_MIN = %d;\n' % K_MIN)
        f.write('constexpr uint32_t NUM_POSSIBLE_TRANSITIONS = %d;\n\n' % ((K_MAX - K_MIN + 1) * (K_MAX + K_MIN + 1)))
        f.write('#endif')

writeHeader()
