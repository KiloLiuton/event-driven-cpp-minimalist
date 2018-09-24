#!/usr/bin/env python
import subprocess

Ks = [5, 10, 20, 40, 80, 160]
n = len(Ks)

for K in Ks:
    N = 10000
    p = 0.0
    subprocess.run(('make N=%d K=%d p=%6.6f' % (N, K, p)).split())
