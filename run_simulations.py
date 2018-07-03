#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor
import subprocess
import time
import sys

executor = ThreadPoolExecutor(max_workers=4)
a = float(sys.argv[1])
for i in range(1, 9):
    N = 500*i
    K = int(N*a)
    executor.submit(subprocess.run, ['./sim-%5.5d-%4.4d-0_0-1_6-graphseed_42' % (N, K), 'trial', 'batch', 'omega'])
