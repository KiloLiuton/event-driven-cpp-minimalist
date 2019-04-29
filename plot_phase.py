#!/usr/bin/env python
"""Plot trial files passed as arguments"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon
from utils import get_header

def getphasetrialdata(filename):
    df = pd.read_csv(filename, skiprows=3, header=None)
    df.columns = list(range(df.shape[1]-1)) + ['time_elapsed']
    return df

def plot_phase(ax, filename, m, M):
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
    title = """N={} K={} p={} gseed={}\na={} iters={} burn={} ic={} seed={}
    stream={}""".format(N, K, p, gs, a, iters, burn, ic, seed, stream)
    ax.set_title(title)
    df = pd.read_csv(filename, skiprows=3, header=None)
    if m is not None and M is not None:
        df = df.iloc[int(m*len(df)):int(M*len(df))]
    df.columns = list(range(df.shape[1]-1)) + ['time_elapsed']
    ax.set_xlabel("Phases")
    ax.set_ylabel("Time [a.u.]")
    print('extent', 0, N, df.time_elapsed.iloc[0], df.time_elapsed.iloc[-1])
    ax.imshow(
        df[df.columns[:-1]].values,
        aspect='auto', origin='lower', cmap='viridis',
        extent=(0, N, df.time_elapsed.iloc[0], df.time_elapsed.iloc[-1])
    )

def main():
    try:
        m, M = sys.argv[-1].split(',')
        m = float(m)
        M = float(M)
        sys.argv = sys.argv[:-1]
    except:
        m,M = None, None
    for filename in sys.argv[1:]:
        fig, ax = plt.subplots(figsize=(10, 8))
        plot_phase(ax, filename, m, M)
        print('Plotted', filename)
    plt.show()

if __name__=="__main__":
    main()
