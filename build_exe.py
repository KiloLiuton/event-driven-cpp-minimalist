#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor
import subprocess
import re
import numpy as np

Ks = [5, 10, 20, 40, 80, 160]
n = len(Ks)
x = np.logspace(-4, -1, n)

for K in Ks:
    N = 10000
    p = 0.0
    subprocess.run(['make', 'N=%d' % N, 'K=%d' % K, 'p=%6.6f' % p])

executor = ThreadPoolExecutor(max_workers=8)
for i in range(n):
    N = 10000
    K = Ks[i]
    p = 0.0
    p = re.sub('\.', '_', '%6.6f'%p)
    cmd = [
        './sim-%5.5d-%4.4d-%s-graphseed_42' % (N, K, p),
        'trial',
        '2.2',
        str(666*i)
    ]
    print(' '.join(cmd))
    subprocess.run(cmd)

