#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from utils import load_batch
from colors import _get_colors

plt.style.use('ggplot')

def get_a_r_psi(folder):
    if not folder.endswith('/'): folder += '/'
    files = sorted([f for f in os.listdir(folder) if f.startswith('batch-N')])
    a = None
    l = []
    r_series = []
    psi_series = []
    for f in files:
        h, data = load_batch(folder + f)
        if a is None:
            a = data.coupling.values
        elif not all(a == data.coupling.values):
            print("coupling values mismatch at file", f)
            sys.exit()
        l.append(1 / h['N'])
        r_series.append(data.r.values)
        psi_series.append(data.psi.values)
    l = np.array(l)
    r_series = np.vstack(r_series).transpose()
    psi_series = np.vstack(psi_series).transpose()
    return l, a, r_series, psi_series

def main():
    folder = sys.argv[1]
    files = [f for f in os.listdir(folder) if f.startswith('batch-N')]
    alpha = 0
    for f in files:
        n = int(f[f.find('N-')+2: f.find('K-')])
        k = int(f[f.find('K-')+2: f.find('p-')])
        alpha += k/n
    alpha /= len(files)
    l, coupling, r_series, psi_series = get_a_r_psi(folder)
    print(l.shape, coupling.shape, r_series.shape, psi_series.shape)
    fig = plt.figure(figsize=(15, 8))

    axr = fig.add_subplot(121)
    axr.set_title('r')
    axr.set_xlabel('$\\lambda$')

    axpsi = fig.add_subplot(122)
    axpsi.set_title('psi')
    axpsi.set_xlabel('$\\lambda$')

    colors = _get_colors(len(r_series))
    for c, a, r, psi in zip(colors, coupling, r_series, psi_series):
        axr.plot(l, r, label=a, color=c)
        axpsi.plot(l, psi, label=f"{a:.3f}", color=c)
    plt.legend(loc=(1.01, 0), title="$a$")
    plt.text(0, 1, f"$\\alpha={alpha:.4f}$", transform=axr.transAxes)
    # plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

