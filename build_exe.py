#!/usr/bin/env python
import numpy as np
import subprocess
import sys

alpha = float(sys.argv[1])
n = 4
x = np.logspace(-4, -1, n)

for i in range(1, n+1):
    N = 3000
    K = int(alpha*N)
    p = x[i-1]
    print('p=%f'%p)
    subprocess.run(['make', 'N=%d'%N, 'K=%d'%K, 'p=%f'%p])
