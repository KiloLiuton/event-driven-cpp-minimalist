#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from colors import _get_colors
import sys
import os


def stacked_plot(ax, X, Y, styles, log=False):
    X = np.array(X)
    Y = np.array(Y)
    n_colors = X.shape[0]
    colors = _get_colors(n_colors)
    for i in range(n_colors):
        if log:
            y = np.log10(Y[i])
        else:
            y = Y[i]
        ax.plot(X[i], y, color=colors[i], *styles[i])


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s ALPHA'.format(sys.argv[0]))
        print('where ALPHA = K/N is the percentage you wish to plot.')
        sys.exit()
    ALPHA = int(100 * float(sys.argv[1]))

    all_batches = sorted(os.listdir('batches'))
    batch_files = []
    for b in all_batches:
        N = int(b[b.find('N-') + 2: b.find('K-')])
        K = int(b[b.find('K-') + 2: b.find('p-')])
        alpha = int(100 * K / N)
        if (alpha == ALPHA):
            batch_files.append(b)
            print(N, K, alpha)

    all_r = []
    all_psi = []
    all_X = []
    all_a = []
    styles_r = []
    styles_psi = []
    styles_X = []
    for f in batch_files:
        N = int(f[f.find('N-') + 2: f.find('K-')])
        batch = np.loadtxt('batches/' + f, skiprows=2, delimiter=',')
        r = batch[:,0]
        r2 = batch[:,1]
        psi = batch[:,2]
        a = batch[:,3]
        X = N * (r2 - r**2)
        sr = ['o']
        spsi = ['o']
        sX = ['o']
        all_r.append(r)
        all_psi.append(psi)
        all_X.append(X)
        all_a.append(a)
        styles_r.append(s)
        styles_psi.append(s)
        styles_X.append(s)
    n_plots = 3
    fig = plt.figure(figsize=(12,10))
    axes = fig.subplots(nrows=n_plots, ncols=1, sharex=True)
    stacked_plot(axes[0], all_a, all_r, styles_r)
    stacked_plot(axes[1], all_a, all_psi, styles_psi)
    stacked_plot(axes[2], all_a, all_X, styles_X, log=True)
    plt.show()
