#!/usr/bin/env python
import pandas as pd
import numpy as np
import os


def get_header(fn, nrows=2):
    with open(fn, "r") as f:
        h = ''
        for _ in range(nrows):
            h = h + f.readline()
    h = h.split()
    d = dict()
    gseed = False
    for x in h:
        if '=' in x:
            key, val = x.split('=')
            key, val = key.strip(), val.strip()
            if key == 'seed' and not gseed:
                key = 'gseed'
                gseed = True
            d[key] = val
    return d


def load_batch(filename):
    ''' Load a batch file and return its contents in a header dict and a
    pandas DataFrame.
    Returns:
    dict(), pd.DataFrame()
    '''
    # TODO: extract header information from header keys in order to improve
    #       mantainability
    # keys = ['N', 'K', 'p', 'seed', 'trials', 'iters', 'burn']
    with open(filename, 'r') as f:
        graph_params = f.readline()
        dynamics_params = f.readline()
    _, N, K, p, seed = graph_params.split(' ')
    _, trials, iters, burn, ic = dynamics_params.split(' ')

    N = int(N.split('=')[1])
    K = int(K.split('=')[1])
    p = float(p.split('=')[1])
    seed = int(seed.split('=')[1])
    trials = int(trials.split('=')[1])
    iters = int(iters.split('=')[1])
    burn = int(burn.split('=')[1])
    ic = ic.split('=')[1].strip()
    headers = dict(
            N=N, K=K, p=p,
            seed=seed,
            trials=trials, iters=iters, burn=burn, initial_condition=ic
            )

    df = pd.read_csv(filename, header=2, delimiter=',')
    return headers, df


def load_trial(filename):
    ''' Load a trial file and return its contents in a header dict and a
    pandas DataFrame.
    Returns:
    dict(), pd.DataFrame()
    '''
    h = get_header(filename)
    h['N'] = int(h['N'])
    h['K'] = int(h['K'])
    h['p'] = float(h['p'])
    h['gseed'] = int(h['gseed'])
    h['coupling'] = float(h['coupling'])
    h['iters'] = int(h['iters'])
    h['burn'] = int(h['burn'])
    h['seed'] = int(h['seed'])
    h['stream'] = int(h['stream'])
    h['initial_condition'] = h['initial_condition'].strip()
    data = pd.read_csv(filename, header=2, delimiter=',')
    return h, data


def makename(N, K, p=0., a0=2., a1=1.35, na=8, v=0):
    """ Create an output file name for a batch simulation given the
    parameters.
    Arguments:
    N -- size of the ring graph
    K -- half the number of neighbors per site in a ring graph

    Keyword arguments:
    p  -- rewiring probability (default 0.0)
    a0 -- minimum coupling stregnth (default 2.0)
    a1 -- maximum coupling strength (default 1.35)
    na -- number os coupling strengths from a0 to a1 (default 8)
    v  -- file version (default 0)
    """
    t = ('batch-N-{0:05d}K-{1:05d}p-{2:6.6f}'
         'a-{3:6.6f}-{4:6.6f}-{5:d}_v{6:d}')
    return t.format(N, K, p, a0, a1, na, v).replace('.', '_') + '.dat'


def plot_batch(axes, data, ylog=True, **plotkwargs):
    ax1, ax2, ax3, ax4 = axes.flatten()
    a = data['coupling']
    # r, w = data['r'], data['omega']
    r = data['r']
    psi = data['psi']
    chi_r, chi_psi = data['chi_r'], data['chi_psi']
    l1 = ax1.plot(a, r, **plotkwargs)
    # ax1.plot(a, w, ls+marker, mfc=mfc, ms=ms)
    l2 = ax2.plot(a, psi, **plotkwargs)
    if ylog:
        l3 = ax3.semilogy(a, chi_r, **plotkwargs)
        l4 = ax4.semilogy(a, chi_psi, **plotkwargs)
    else:
        l3 = ax3.plot(a, chi_r, **plotkwargs)
        l4 = ax4.plot(a, chi_psi, **plotkwargs)
    return l1, l2, l3, l4


def getFilesWithAlpha(folder, alpha):
    files = [f for f in os.listdir(folder) if 'batch-' in f]
    nlist = []
    for f in files:
        i = f.find("N-") + 2
        N = int(f[i: f.find("K-")])
        nlist = nlist + [N] if not N in nlist else nlist
    nlist.sort(reverse=True)
    result = []
    for N in nlist:
        K = int(round(N*alpha))
        result += [f for f in files if 'N-%05d'%N in f and 'K-%05d'%K in f]
    return result
    
