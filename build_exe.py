#!/usr/bin/env python
import subprocess
import sys

a = float(sys.argv[1])
p = 0.0

for i in range(1, 9):
    N = 500*i
    K = int(a*N)
    subprocess.run(['make', 'N=%d'%N, 'K=%d'%K, 'p=%.1f'%p])
