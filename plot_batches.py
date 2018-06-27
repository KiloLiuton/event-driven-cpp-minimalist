#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

for filename in sys.argv[1:]:
    N = filename[filename.find('N-')+2:filename.find('K-')]
    K = filename[filename.find('K-')+2:filename.find('p-')]
    title = 'N=' + N + ' K=' + K

    data = np.loadtxt(filename, skiprows=2, delimiter=',')

    r = data[:,0]
    r2 = data[:,1]
    psi = data[:,2]
    a = data[:,3]

    print('foooo', type(r2), type(r**2))

    X = int(N)*(r2 - r**2)

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax.set_title(title)

    ax.set_xlabel('a')
    ax.set_ylabel('r')
    ax2.set_xlabel('a')
    ax2.set_ylabel('X')
    ax3.set_xlabel('a')
    ax3.set_ylabel('psi')


    ax.plot(a, r, 'ro')
    ax2.plot(a, X, 'go')
    ax3.plot(a, psi, 'bo')

plt.show()
