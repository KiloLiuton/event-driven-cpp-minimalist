import sys
import random

try:
    N, K, p = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
except:
    sys.exit('Usage:\n./generate_header.py N K p\nN=size, K=forward neighbors, p=rewire probability')
if 2*K > N:
    sys.exit('N must be greater than 2*K. Got N=%d, K=%d'%(N, K))

def printHeader():
    with open('header.h', 'w') as f:
        # TODO: print the actual content
        f.write('foo!')

# REGULAR RING TOPOLOGY
NEIGHBOR_LIST = [0]*(N*2*K)
t = 0
for i in range(N):
    for j in range(i - K, i + K + 1):
        if j == i: continue
        NEIGHBOR_LIST[t] = j % N
        t += 1

INDEXES = [2*K*i for i in range(N)]
NUMBER_OF_NEIGHBORS = [2*K]*N
K_MAX = max(NUMBER_OF_NEIGHBORS)
K_MIN = min(NUMBER_OF_NEIGHBORS)
NUM_POSSIBLE_TRANSITIONS = (K_MAX - K_MIN + 1) * (K_MAX + K_MIN + 1)

print(INDEXES, len(INDEXES))
print(NEIGHBOR_LIST, len(NEIGHBOR_LIST))
print(NUMBER_OF_NEIGHBORS, len(NUMBER_OF_NEIGHBORS))

class ListDict(object):
    def __init__(self, l = []):
        self.item_to_position = {}
        self.items = l

    def add_item(self, item):
        if item in self.item_to_position:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove_item(self, item):
        position = self.item_to_position.pop(item)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position

    def choose_random_item(self):
        return random.choice(self.items)

# REWIRE
if p == 0:
    printHeader()
else:
    all = set(range(N))
    for j in range(K, 2*K):
        for i in range(N):
            r = random.random()
            if r > p:
                continue
            i1, i2 = INDEXES[i], INDEXES[i + NUMBER_OF_NEIGHBORS[i]]
            choices = all - set(NEIGHBOR_LIST[i1:i2] + [i])
    K_MAX = max(NUMBER_OF_NEIGHBORS)
    K_MIN = min(NUMBER_OF_NEIGHBORS)
    NUM_POSSIBLE_TRANSITIONS = (K_MAX - K_MIN + 1) * (K_MAX + K_MIN + 1)
    printHeader();

print('DONE!')
