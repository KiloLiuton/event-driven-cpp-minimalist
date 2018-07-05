#!/usr/bin/env python
import numpy as np
import subprocess
import sys

alpha = float(sys.argv[1])
n = 5
x = np.logspace(-4, -1, n)

for i in range(1, n+1):
    N = 500*i
    K = int(alpha*N)
    p = 0.0
    print('p=%f'%p)
    subprocess.run(['make', 'N=%d'%N, 'K=%d'%K, 'p=%6.6f'%p])
