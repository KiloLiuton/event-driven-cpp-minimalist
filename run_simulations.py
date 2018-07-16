#!/usr/bin/env python
from concurrent.futures import ThreadPoolExecutor
import subprocess
import time
import sys
import re
import numpy as np

executor = ThreadPoolExecutor(max_workers=8)
n = 10
for i in range(1, n+1):
    N = 2000
    K = 680
    p = 0.0
    p = re.sub('\.', '_', '%6.6f'%p)
    cmd = ['./sim-%5.5d-%4.4d-%s-graphseed_42' % (N, K, p), 'trial', '2.2', str(666*i)]
    print(' '.join(cmd))
    subprocess.run(cmd)
    #executor.submit(subprocess.run, cmd)
