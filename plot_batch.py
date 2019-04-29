#!/usr/bin/env python
"""Plot batch files passed as arguments"""
import sys

from utils import load_batch, plot_batch
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from colors import _get_colors


def ls(i):
    styles = ['-', ':', '-.', '--']
    return styles[i % len(styles)]


def ms(i):
    styles = ['s', 'v', 'o', 'P']
    return styles[i % len(styles)]


fs = 12
fig = plt.figure(figsize=(fs, fs*9/16))
axes = fig.subplots(2, 2, sharex=True)
ax1, ax2, ax3, ax4 = axes.flatten()

# ax1.legend(loc='lower right', fontsize=fs)
ax1.set_ylabel('$r$', fontsize=fs*1.4)
ax2.set_ylabel('$\\psi$', fontsize=fs*1.4)
ax3.set_ylabel('$\\chi_r$', fontsize=fs*1.4)
ax4.set_ylabel('$\\chi_{\\psi}$', fontsize=fs*1.4)

ax1.text(0.02, 0.94, 'a)', transform=ax1.transAxes, fontsize=fs)
ax2.text(0.02, 0.94, 'b)', transform=ax2.transAxes, fontsize=fs)
ax3.text(0.02, 0.94, 'c)', transform=ax3.transAxes, fontsize=fs)
ax4.text(0.02, 0.94, 'd)', transform=ax4.transAxes, fontsize=fs)

ax3.set_xlabel('Coupling Strength', fontsize=fs)
ax4.set_xlabel('Coupling Strength', fontsize=fs)
for ax in axes.flatten():
    ax.tick_params(labelsize=fs)

handles = []
labels = []
ic = ''
for i, arg in enumerate(sys.argv):
    if arg.startswith("--step="):
        step = int(sys.argv.pop(i).split('=')[1])
        break
    else:
        step = 1
colors = _get_colors(int((len(sys.argv) - 1)/step))
for i, filename in enumerate(sys.argv[1:][::step]):
    header, data = load_batch(filename)

    if not ic:
        ic = header['initial_condition']
    elif ic != header['initial_condition']:
        ic = 'mixed'

    N, K = header['N'], header['K']
    data['chi_r'] = data['chi_r']*N
    data['chi_psi'] = data['chi_psi']*N
    print(header)
    plotkwargs = dict(color=colors[i], marker=ms(i), mfc='w', ms=fs*0.5)
    plot_batch(axes, data, ylog=True, **plotkwargs)
    handles += [Line2D([], [], **plotkwargs)]
    alpha = K / N
    if header['p']:
        title = 'N  $\\alpha$  p'
        tmp = '{:04d}  {:f}  {:f}'.format(header['N'], alpha, header['p'])
    else:
        title = 'N             $\\alpha$'
        tmp = '{:04d}  {:f}'.format(header['N'], alpha)
    labels += [tmp]
plt.figlegend(
        handles,
        labels,
        fontsize=fs,
        bbox_to_anchor=(1.0, 1.0),
        bbox_transform=fig.transFigure,
        title=title,
        title_fontsize=fs
        )

fig.suptitle('Initial Condition: {}'.format(ic))
fig.tight_layout(rect=[0, 0, 0.84, 0.95])
plt.show()
