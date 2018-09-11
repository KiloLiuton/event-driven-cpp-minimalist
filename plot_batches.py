#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import colors
import sys

for filename in sys.argv[1:]:
    N = filename[filename.find('N-')+2:filename.find('K-')]
    K = filename[filename.find('K-')+2:filename.find('p-')]
    p = filename[filename.find('p-')+2:filename.find('a-')]
    title = 'N=' + N + ' K=' + K + ' p=' + p

    data = pd.read_csv(filename, header=2, delimiter=',')
    print(data.head())

    r = data['r']
    # coupling,r,r2,psi,psi2,chi_r,chi_psi,omega,used_seed
    a = data['coupling'].values
    r2 = data['r2'].values
    psi = data['psi'].values
    psi2 = data['psi2'].values
    chi_r = data['chi_r'].values
    chi_psi = data['chi_psi'].values
    omega = 2 * np.pi * data['omega'].values
    batch_seeds = data['used_seed'].values

    chi_r = int(N)*(r2 - r**2)

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax.set_title(title)

    ax.set_xlabel('a')
    ax.set_ylabel('r')
    ax2.set_xlabel('a')
    ax2.set_ylabel('$\\chi$')
    ax3.set_xlabel('a')
    ax3.set_ylabel('psi')

    clist = colors._get_colors(4)
    ax.plot(a, r, 'o', color=clist[0])
    ax.plot(a, omega, 'o', color=clist[1])
    ax2.semilogy(a, chi_r, 'o', color=clist[2])
    ax3.plot(a, psi, 'o', color=clist[3])

plt.show()
