#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor
import subprocess
import time
import sys
import re
import numpy as np

executor = ThreadPoolExecutor(max_workers=4)
#alpha = float(sys.argv[1])
a = float(sys.argv[1])
x = np.logspace(-4, 0, 8)
for i in range(1, 9):
    N = 3000
    K = 30
    p = re.sub('\.', '_', '%6.6f'%x[i-1])
    executor.submit(subprocess.run, ['./sim-%5.5d-%4.4d-%s-graphseed_42' % (N, K, p), 'trial', '%f'%a])
