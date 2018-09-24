import subprocess
import numpy as np
import os
np.set_printoptions(precision=3)

num_points = 8
max_N = 100
min_N = 30
num_couplings = 8
max_coupling = 1.8
min_coupling = 1.35
alpha = 0.34
p = 0.0

N_list = np.flipud(1/np.linspace(1/max_N, 1/min_N, num_points)).tolist()
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
