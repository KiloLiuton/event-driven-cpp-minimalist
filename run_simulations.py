#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor
import subprocess
import time
import sys
import re
import numpy as np

executor = ThreadPoolExecutor(max_workers=4)
alpha = float(sys.argv[1])
n = 4
print('foo')
for i in range(1, n+1):
    N = 500 * i
    K = int(N * alpha)
    p = 0.0
    p = re.sub('\.', '_', '%6.6f'%p)
    print(' '.join(['./sim-%5.5d-%4.4d-%s-graphseed_42' % (N, K, p), 'trial', '2.0', 'batch', 'omega']))
    executor.submit(subprocess.run, ['./sim-%5.5d-%4.4d-%s-graphseed_42' % (N, K, p), 'trial', '2.0', 'batch', 'omega'])
