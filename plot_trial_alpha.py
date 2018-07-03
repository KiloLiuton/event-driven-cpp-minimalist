#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from colors import _get_colors
import sys
import os


def stacked_plot(ax, X, Y, log=False, styling = ['o']):
    X = np.array(X)
    Y = np.array(Y)
    n_colors = X.shape[0]
    colors = _get_colors(n_colors)
    for i in range(n_colors):
        if log:
            y = np.log10(Y[i])
        else:
            y = Y[i]
        ax.plot(X[i], y, color=colors[i], **styling)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: {} ALPHA FOLDER'.format(sys.argv[0]))
        print('where ALPHA = K/N is the percentage you wish to plot.')
        print('where FOLDER is the path to the trial data.')
        sys.exit()
    ALPHA = int(100 * float(sys.argv[1]))
    FOLDER = sys.argv[2]
    if not FOLDER.endswith('/'): FOLDER = FOLDER + '/'

    all_trials = sorted([f for f in os.listdir(FOLDER) if not f.startswith('batch')])
    print('BARRRR', all_trials)
    trial_files = []
    for t in all_trials:
        N = int(t[t.find('N-') + 2: t.find('K-')])
        K = int(t[t.find('K-') + 2: t.find('p-')])
        alpha = int(100 * K / N)
        if (alpha == ALPHA):
            trial_files.append(t)
            print(N, K, alpha)

    all_r = []
    all_psi = []
    all_X = []
    all_a = []
    for t in trial_files:
        N = int(t[t.find('N-') + 2: t.find('K-')])
        batch = np.loadtxt(FOLDER + t, skiprows=2, delimiter=',')
        r = batch[:,0]
        r2 = batch[:,1]
        psi = batch[:,2]
        a = batch[:,3]
        X = N * (r2 - r**2)
        all_r.append(r)
        all_psi.append(psi)
        all_X.append(X)
        all_a.append(a)
    n_plots = 3
    fig = plt.figure(figsize=(12,10))
    axes = fig.subplots(nrows=n_plots, ncols=1, sharex=True)
    stacked_plot(axes[0], all_a, all_r, styling=['-'])
    stacked_plot(axes[1], all_a, all_psi, styling=['-'])
    stacked_plot(axes[2], all_a, all_X, log=True, styling=['-'])
    plt.show()
