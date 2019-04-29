#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from utils import load_trial, get_header

def plot_phase(ax, filename):
    h = get_header(filename)
    N = int(h['N'])
    K = int(h['K'])
    p = float(h['p'])
    gs = int(h['gseed'])
    a = float(h['coupling'])
    iters = int(h['iters'])
    burn = int(h['burn'])
    ic = h['initial_condition']
    seed = int(h['seed'])
    stream = int(h['stream'])
    title = """N={}  K={}  p={}  gseed={}  a={}  ic={}\niters={}  burn={}  seed={}\
  stream={}""".format(N, K, p, gs, a, ic, iters, burn, seed, stream)
    # ax.set_title(title)
    df = pd.read_csv(filename, skiprows=3, header=None)
    df.columns = list(range(df.shape[1]-1)) + ['time_elapsed']
    ax.set_xlabel("Phases", fontsize=14)
    ax.set_ylabel("Time [a.u.]", fontsize=14)
    # CROP PART OF THE TIME AXIS
    nrows = df.shape[0]
    s, e = 0, -1
    M = df.iloc[s:e]
    ax.imshow(
            M.values[:,0:-1],
            aspect='auto', origin='lower', cmap='viridis',
            extent=(0, N, M.time_elapsed.iloc[0], M.time_elapsed.iloc[-1])
            )
    return title

def plotpops(ax, p0, p1, t):
    p0 = [(t[0], 0)] + list(zip(t, p0 + p1)) + [(t[-1], 0)]
    p1 = [(t[0], 0)] + list(zip(t, p1)) + [(t[-1], 0)]
    poly1 = Polygon(p0, facecolor='magenta', edgecolor='k')
    poly2 = Polygon(p1, facecolor='lawngreen', edgecolor='k')
    ax.set_facecolor('yellow')
    ax.add_patch(poly1)
    ax.add_patch(poly2)
    ax.plot([], [], 's', color='magenta', label="$N_0$")
    ax.plot([], [], 's', color='lawngreen', label="$N_1$")
    ax.plot([], [], 's', color='yellow', label="$N_2$")
    ax.legend(facecolor='w', framealpha=.75, fontsize=14, markerscale=2)

def usage():
    sys.exit("Usage: Choose one of the two options\nPlot an identical" +
             "corresponding file: ./plot_pair [trialfile OR" +
             "phasefile]\nExplicitly choose both files: ./plot_pair [trialfile]" +
             "[phasefile]\n")

tfile = pfile = ''
if len(sys.argv) == 2:
    tfile = sys.argv[1]
    if tfile.startswith("phase"):
        pfile = tfile
        tfile = tfile[5:]
    else:
        pfile = "phase" + tfile
elif len(sys.argv) == 3:
    tfile = sys.argv[1]
    pfile = sys.argv[2]
    if tfile.startswith("phase"):
        tfile, pfile = pfile, tfile
else:
    usage()

try:
    h, data = load_trial(tfile)
except:
    print("could not oprn file", tfile)
    usage()
N = h['N']
K = h['K']
p = h['p']
gs = h['gseed']
a = h['coupling']
iters = h['iters']
burn = h['burn']
ic = h['initial_condition']
seed = h['seed']
stream = h['stream']
xmin, xmax = 0.00, 1.00
L = data.shape[0]
r = np.sqrt(data['r'].values)[int(xmin*L): int(xmax*L)]
psi = data['psi'].values[int(xmin*L): int(xmax*L)]
p0 = data['pop0'].values[int(xmin*L): int(xmax*L)] / N
p1 = data['pop1'].values[int(xmin*L): int(xmax*L)] / N
t = data['time_elapsed'].values[int(xmin*L): int(xmax*L)]
dt = data['dt'].values[int(xmin*L): int(xmax*L)]
filesize = os.path.getsize(tfile)
if filesize > 300e6:
    red = int(filesize / 300e6)
    r = r[::red]
    psi = psi[::red]
    p0 = p0[::red]
    p1 = p1[::red]
    t = t[::red]
    dt = dt[::red]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(121)
captionleft = plot_phase(ax, pfile)
ax.tick_params(labelsize=14)

ax1 = fig.add_subplot(322)
captionright = """N={}  K={}  p={}  gseed={}  a={}  ic={}\niters={}  burn={}  seed={}\
  stream={}""".format(N, K, p, gs, a, ic, iters, burn, seed, stream)
if captionleft == captionright:
    fig.text(.4, .05, captionleft, fontsize=14)
else:
    fig.text(.6, .05, captionleft, fontsize=14)
    fig.text(.1, .05, captionright, fontsize=14)
ax1.set_xlim(t[0], t[-1])
ax1.set_ylim(0, 1)
ax1.tick_params(labelsize=14)
ax1.set_ylabel('Populations', fontsize=14)
plotpops(ax1, p0, p1, t)

ax2 = fig.add_subplot(324, sharex=ax1)
ax2.set_xlim(t[0], t[-1])
ax2.set_ylim(0, 1)
ax2.tick_params(labelsize=14)
ax2.set_ylabel('$r$', fontsize=14)
ax2.set_ylim(0, max(r)*1.15)
ax2.plot(t, r, 'b-')

ax3 = fig.add_subplot(326, sharex=ax1)
ax3.set_xlim(t[0], t[-1])
ax3.set_ylim(0, max(psi)*1.15)
ax3.tick_params(labelsize=14)
ax3.set_xlabel('Time [a.u.]', fontsize=14)
ax3.set_ylabel('$\\psi$', fontsize=14)
ax3.plot(t, psi, 'r-')

if (burn < xmax*L) and (burn > xmin*L):
    tburn = t[int((burn - int(xmin*L))/red)]
    ax.axhline(tburn, color='red', lw=1.5, linestyle='--')
    ax1.axvline(tburn, color='red', lw=1.5, linestyle='--')
    ax2.axvline(tburn, color='red', lw=1.5, linestyle='--')
    ax3.axvline(tburn, color='red', lw=1.5, linestyle='--')

plt.tight_layout()
plt.tight_layout(rect=[0, .085, 1, 1])
print('Plotted:', pfile, tfile)

plt.show()
