#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import colors
import sys

for filename in sys.argv[1:]:
    N = filename[filename.find('N-')+2:filename.find('K-')]
    K = filename[filename.find('K-')+2:filename.find('p-')]
    title = 'N=' + N + ' K=' + K

    data = np.loadtxt(filename, skiprows=2, delimiter=',')

    r = data[:, 0]
    r2 = data[:, 1]
    psi = data[:, 2]
    omega = 2 * np.pi * data[:, 3]
    omega[:8] = 0
    a = data[:, 4]

    X = int(N)*(r2 - r**2)

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax.set_title(title)

    ax.set_xlabel('a')
    ax.set_ylabel('r')
    ax2.set_xlabel('a')
    ax2.set_ylabel('X')
    ax3.set_xlabel('a')
    ax3.set_ylabel('psi')

    clist = colors._get_colors(4)
    ax.plot(a, r, 'o', color=clist[0])
    ax.plot(a, omega, 'o', color=clist[1])
    ax2.semilogy(a, X, 'o', color=clist[2])
    ax3.plot(a, psi, 'o', color=clist[3])

plt.show()
