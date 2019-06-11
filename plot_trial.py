#!/usr/bin/env python
"""Plot trial files passed as arguments"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from utils import load_trial

for filename in sys.argv[1:]:
    h, data = load_trial(filename)

    N = h['N']
    K = h['K']
    p = h['p']
    a = h['coupling']
    iters = h['iters']
    burn = h['burn']
    ic = h['initial_condition']
    seed = h['seed']
    stream = h['stream']
    title = 'N={} K={} p={} a={} iters={} burn={}\nic={} seed={} stream={}'.format(
            N, K, p, a, iters, burn, ic, seed, stream)
    fs = 16
    xmin = 0.00
    xmax = 1.00
    numlines = data.shape[0]
    r = np.sqrt(data['r'].values)[int(xmin*numlines): int(xmax*numlines)]
    psi = data['psi'].values[int(xmin*numlines): int(xmax*numlines)]
    p0 = data['pop0'].values[int(xmin*numlines): int(xmax*numlines)] / N
    p1 = data['pop1'].values[int(xmin*numlines): int(xmax*numlines)] / N
    t = data['time_elapsed'].values[int(xmin*numlines): int(xmax*numlines)]
    dt = data['dt'].values[int(xmin*numlines): int(xmax*numlines)]
    # if filesize is > 300MB, discard one line every `red` lines
    filesize = os.path.getsize(filename)
    red = 1
    if filesize > 300e6:
        red = int(filesize / 300e6)
        r = r[::red]
        psi = psi[::red]
        p0 = p0[::red]
        p1 = p1[::red]
        t = t[::red]
        dt = dt[::red]
    fig = plt.figure(figsize=(fs, fs*9/16))
    ax1 = fig.add_subplot(311, facecolor='yellow')
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax1)
    ax1.set_title(title, fontsize=fs*1.6)
    ax2.set_ylim(0, 1)
    ax3.set_ylim(0, 1)

    ax1.tick_params(labelsize=fs)
    ax2.tick_params(labelsize=fs)
    ax3.tick_params(labelsize=fs)

    tmp = int(numlines*(burn/iters - xmin)/red)
    if tmp < 0: tmp = 0
    ax1.axvline(t[tmp], color='red', lw=2)
    ax2.axvline(t[tmp], color='red', lw=2)

    ax1.set_xlim(t[0], t[-1])
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Populations', fontsize=fs*1.4)
    ax2.set_ylabel('$r$', fontsize=fs*1.4)
    ax2.set_xlim(t[0], t[-1])
    ax2.set_ylim(0, 1)
    ax3.set_xlabel('Time [a.u.]', fontsize=fs*1.4)
    ax3.set_ylabel('$\\psi$', fontsize=fs*1.4)
    ax3.set_xlim(t[0], t[-1])
    box = dict(facecolor='gray', alpha=.5)
    ax3.text(.9, .8, '$<\\psi>$=%f\n$<r>=$%f'%(np.mean(psi[tmp:]), np.mean(r[tmp:])), bbox=box, transform=ax3.transAxes, fontsize=fs)


    p0 = [(t[0], 0)] + list(zip(t, p0 + p1)) + [(t[-1], 0)]
    p1 = [(t[0], 0)] + list(zip(t, p1)) + [(t[-1], 0)]
    poly1 = Polygon(p0, facecolor='magenta', edgecolor='k')
    poly2 = Polygon(p1, facecolor='lawngreen', edgecolor='k')

    ax1.add_patch(poly1)
    ax1.add_patch(poly2)
    ax1.plot([], [], 's', color='magenta', label="$N_0$")
    ax1.plot([], [], 's', color='lawngreen', label="$N_1$")
    ax1.plot([], [], 's', color='yellow', label="$N_2$")
    ax1.legend(facecolor='w', framealpha=.75, fontsize=fs/1.4, markerscale=fs/8)
    ax2.plot(t, r, 'b-')
    ax3.plot(t, psi, 'r-')
    plt.tight_layout()
    print('Plotted', filename)

plt.show()
