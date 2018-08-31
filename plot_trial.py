#!/usr/bin/env python
"""Plot a trial files passed as arguments"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon

for filename in sys.argv[1:]:
    N = filename[filename.find('N-')+2:filename.find('K-')]
    K = filename[filename.find('K-')+2:filename.find('p-')]
    title = 'N=' + N + ' K=' + K

    datap = pd.read_csv(filename, header=1, delimiter=',')
    print(datap.head())

    with open(filename, 'r') as f:
        header = f.readline()
    N, K, p, a, BURN, ITERS = header.split(',')
    N = int(N.split('=')[1])
    K = int(K.split('=')[1])
    p = float(p.split('=')[1])
    a = float(a.split('=')[1])
    BURN = int(BURN.split('=')[1])
    ITERS = int(ITERS.split('=')[1])

    xmin = 0.25
    xmax = 1
    L = datap.shape[0]
    r = np.sqrt(datap['r**2'].values)[int(xmin*L): int(xmax*L)]
    psi = datap['psi'].values[int(xmin*L): int(xmax*L)]
    n0 = datap['N0'].values[int(xmin*L): int(xmax*L)]
    n1 = datap['N1'].values[int(xmin*L): int(xmax*L)]
    t = datap['time'].values[int(xmin*L): int(xmax*L)]

    fig = plt.figure()
    ax = fig.add_subplot(211, facecolor='yellow')
    ax2 = fig.add_subplot(212)
    ax.set_title(title)
    ax2.set_ylim(0, 1)

    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, N)
    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Order parameter r')

    # ax.axvline(t[BURN], color='red', lw=2)
    # ax2.axvline(t[BURN], color='red', lw=2)

    n0 = [(t[0], 0)] + list(zip(t, n0 + n1)) + [(t[-1], 0)]
    n1 = [(t[0], 0)] + list(zip(t, n1)) + [(t[-1], 0)]
    poly1 = Polygon(n0, facecolor='magenta', edgecolor='k')
    poly2 = Polygon(n1, facecolor='lawngreen', edgecolor='k')

    ax.add_patch(poly1)
    ax.add_patch(poly2)
    ax2.plot(t, r, 'b-')

# plt.savefig('foooo.png')
plt.show()

