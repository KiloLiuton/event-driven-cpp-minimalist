#!/usr/bin/env python
import subprocess
import numpy as np
import sys
import os
np.set_printoptions(precision=3)


def usage():
    print('Usage: python %s args' % sys.argv[0])
    print('--num_points   [int] number of system sizes')
    print('--max_N        [int] maximum syste size')
    print('--min_N        [int] minumum syste size')
    print('--alpha        [float] K/N ratio')
    print('--max_coupling [float] maximum coupling strength')
    print('--min_coupling [float] minimum coupling strength')
    print('--num_coupling [int] number of couplings per system size')
    print('--prob         [float] rewiring probability')
    sys.exit()


num_points = 8
max_N = 500
min_N = 100
num_couplings = 8
max_coupling = 1.8
min_coupling = 1.35
alpha = 0.12
p = 0.000000

# parse command line arguments
for i, arg in enumerate(sys.argv):
    if arg.startswith('-h') or arg.startswith('--help'):
        usage()
    try:
        if arg.startswith('--num_points'):
            num_points = int(sys.argv[i+1])
        if arg.startswith('--max_N'):
            max_N = int(sys.argv[i+1])
        if arg.startswith('--min_N'):
            min_N = int(sys.argv[i+1])
        if arg.startswith('--alpha'):
            alpha = float(sys.argv[i+1])
        if arg.startswith('--max_coupling'):
            max_coupling = float(sys.argv[i+1])
        if arg.startswith('--min_coupling'):
            min_coupling = float(sys.argv[i+1])
        if arg.startswith('--num_couplings'):
            num_couplings = int(sys.argv[i+1])
        if arg.startswith('--prob'):
            p = float(sys.argv[i+1])
    except(ValueError):
        print("Some arguments were not recognized!")
        usage()

N_list = np.flipud(1/np.linspace(1/max_N, 1/min_N, num_points)).astype(int).tolist()
K_list = [int(alpha * n) for n in N_list]
couplings = np.linspace(1.35, 1.8, num_couplings).tolist()

# compile files
print('Compiling files...')
for n, k in zip(N_list, K_list):
    command = 'make N=%d K=%d p=%6.6f' % (n, k, p)
    name = ('sim-%5.5d-%4.4d-%.6f-gseed_42' % (n, k, p)).replace('.', '_')
    if name not in os.listdir():
        subprocess.run(command.split())
    else:
        print('Executable %s already present!' % name)
print('Compiling Done!')

# execute files
print('Executing files...')
params = [
        '-b',
        '-bs', '%.6f' % min_coupling,
        '-be', '%.6f' % max_coupling,
        '-bn', '%d' % num_couplings
    ]
for n, k in zip(N_list, K_list):
    name = ('sim-%5.5d-%4.4d-%.6f-gseed_42' % (n, k, p)).replace('.', '_')
    command = ['./' + name]
    print(' '.join(command + params))
    subprocess.run(command + params)
print('Execution done, Congratulations!')
