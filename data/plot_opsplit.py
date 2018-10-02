#!/usr/bin/env python
''' Concatenate all chi curves experiments to for a op_split graph
'''

import matplotlib.pyplot as plt
import itertools as it
import pandas as pd
import numpy as np
import colors
import os

folders = sorted(
    [f for f in os.listdir() if os.path.isdir(f) and not f == '__pycache__'])
for i, f in enumerate(folders):
    print(i+1, f)
answer = int(input('Which folder do you want to open? ')) - 1
folder = folders[answer] + '/'

files = sorted(
        [f for f in os.listdir(folder) if f.startswith('batch-')])

df = pd.DataFrame(columns=['N', 'K', 'p', 'alpha', 'graph_seed',
                           'trials', 'iters', 'burn'])
df2 = pd.DataFrame()
couplings = []
curves_r = []
curves_psi = []
curves_chi_r = []
curves_chi_psi = []
for f in files:
    f = folder + f
    with open(f, 'r') as ff:
        header1 = ff.readline().strip().split()[1:]
        header2 = ff.readline().strip().split()[1:]
    N, K, p, gseed = [h.split('=')[1] for h in header1]
    trials, iters, burn = [h.split('=')[1] for h in header2]
    N, K, trials, iters, burn = [int(x) for x in [N, K, trials, iters, burn]]
    idx = df.shape[0]
    df.loc[idx] = [N, K, p, K/N, gseed, trials, iters, burn]
    data = pd.read_csv(f, sep=',', header=2)
    if not couplings:
        couplings = data.coupling.tolist()
    curves_r.append(data.r.tolist())
    curves_psi.append(data.psi.tolist())
    curves_chi_r.append(data.chi_r.tolist())
    curves_chi_psi.append(data.chi_psi.tolist())

curves_r = np.array(curves_r).transpose().tolist()
curves_psi = np.array(curves_psi).transpose().tolist()
curves_chi_r = np.array(curves_chi_r).transpose().tolist()
curves_chi_psi = np.array(curves_chi_psi).transpose().tolist()

for a, c1, c2, c3, c4 in zip(
        couplings, curves_r, curves_psi,
        curves_chi_r, curves_chi_psi):
    df2.insert(0, 'r-%6.6f' % a, c1)
    df2.insert(0, 'psi-%6.6f' % a, c2)
    df2.insert(0, 'chi_r-%6.6f' % a, c3)
    df2.insert(0, 'chi_psi-%6.6f' % a, c4)
df2 = df2.reindex(sorted(df2.columns), axis=1)
df = df.join(df2)
alpha = (df.K / df.N).mean()


def plt_series(df, prefix, ylabel, title=''):
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, facecolor=(0.9, 0.9, 0.9,))
    cl = it.cycle(colors._get_colors(8))  # colors list
    mk = colors.markerStyles()            # marker styles list
    maxy = -float('infinity')
    miny = float('infinity')
    for column in reversed(df.columns):
        if column.startswith(prefix):
            color = next(cl)
            ax.plot(1/df.N, np.log10(df[column]),
                    lw=2.3, color=color,
                    marker=mk.next_thick(), ms=16,
                    alpha=0.4, mfc=color, label=column.split('-')[1][:5])
            maxy = max(max(df[column]), maxy)
            miny = min(min(df[column]), miny)
    bbox_props = dict(
            boxstyle="round", fc="w", ec="0.0", alpha=0.55)
    ax.text(0.82, 0.12, "Coupling\nStrength",
            ha="center", va="center", size=20,
            bbox=bbox_props, transform=ax.transAxes)
    ax.text(0.5, 0.12,
            "Parameters:\n$\\alpha=%.2f$\n$p=%.5f$" % (alpha, float(p)),
            ha="center", va="center", size=20,
            bbox=bbox_props, transform=ax.transAxes)
    plt.legend(loc=4, fontsize=14)
    yticks = np.linspace(miny, maxy, 10)
    ax.grid(color='w', which='both')
    e = 0.10            # spacing on top and bottom values of the y-axis
    y0, y1 = np.log10([miny, maxy])
    ax.set_ylim(y0 - e*(y1-y0), y1 + e*(y1-y0))
    ax.set_yticks(np.log10(yticks))
    ax.set_yticklabels(['%.3f' % y for y in yticks], fontsize=15)
    ax.set_xticks(1/df.N)
    ax.set_xticklabels([str(int(n)) for n in df.N], fontsize=15)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_xlabel('System Size $N$', fontsize=20)
    ax.set_title(title, fontsize=18)


plt_series(df, 'r-', 'Order Parameter $r$')
plt_series(df, 'psi-', 'Order Parameter $\\psi$')
plt.show()
# plt.savefig('foo.png')
