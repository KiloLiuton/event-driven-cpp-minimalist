#!/usr/bin/env python
"""Plot a trial files passed as arguments"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon

for filename in sys.argv[1:]:
    datap = pd.read_csv(filename, header=2, delimiter=',')
    print(datap.head())

    with open(filename, 'r') as f:
        graph_params = f.readline()
        dynamics_params = f.readline()
    _, N, K, p, seed = graph_params.split(' ')
    _, a, iters, burn, _, _ = dynamics_params.split(' ')
    title = N + ' ' + K + ' ' + p + ' ' + a
    N = int(N.split('=')[1])
    K = int(K.split('=')[1])
    p = float(p.split('=')[1])
    a = float(a.split('=')[1])
    iters = int(iters.split('=')[1])
    burn = int(burn.split('=')[1])

    xmin = 0.00
    xmax = 1.00
    L = datap.shape[0]
    r = np.sqrt(datap['r'].values)[int(xmin*L): int(xmax*L)]
    psi = datap['psi'].values[int(xmin*L): int(xmax*L)]
    n0 = datap['pop0'].values[int(xmin*L): int(xmax*L)]
    n1 = datap['pop1'].values[int(xmin*L): int(xmax*L)]
    t = datap['time_elapsed'].values[int(xmin*L): int(xmax*L)]
    dt = datap['dt'].values[int(xmin*L): int(xmax*L)]

    fig = plt.figure()
    ax = fig.add_subplot(311, facecolor='yellow')
    ax2 = fig.add_subplot(312, sharex=ax)
    ax3 = fig.add_subplot(313, sharex=ax)
    ax.set_title(title)
    ylim = max(max(psi), max(r))
    ax2.set_ylim(0, ylim)
    ax3.set_ylim(0, ylim)

    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, N)
    ax.set_ylabel('Populations')
    ax2.set_ylabel('$r$')
    ax2.set_xlim(t[0], t[-1])
    ax3.set_xlabel('Time')
    ax3.set_ylabel('$\\psi$')
    ax3.set_xlim(t[0], t[-1])

    if (burn < xmax*L) and (burn > xmin*L):
        ax.axvline(t[burn - int(xmin*L)], color='red', lw=2)
        ax2.axvline(t[burn - int(xmin*L)], color='red', lw=2)

    n0 = [(t[0], 0)] + list(zip(t, n0 + n1)) + [(t[-1], 0)]
    n1 = [(t[0], 0)] + list(zip(t, n1)) + [(t[-1], 0)]
    poly1 = Polygon(n0, facecolor='magenta', edgecolor='k')
    poly2 = Polygon(n1, facecolor='lawngreen', edgecolor='k')

    ax.add_patch(poly1)
    ax.add_patch(poly2)
    ax2.plot(t, r, 'b-')
    ax3.plot(t, psi, 'r-')

plt.tight_layout()
# plt.savefig('foooo.png')
plt.show()

