#!/usr/bin/env python
import numpy as np
import subprocess

alphas = np.linspace(0.01, 0.45, 8)
n = 4
x = np.logspace(-4, -1, n)

for alpha in alphas:
    N = 1001
    K = int(alpha*N)
    p = 0.0
    subprocess.run(['make', 'N=%d' % N, 'K=%d' % K, 'p=%6.6f' % p])

