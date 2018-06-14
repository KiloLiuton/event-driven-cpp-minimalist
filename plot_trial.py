#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

for filename in sys.argv[1:]:
    N = filename[filename.find('N-')+2:filename.find('K-')]
    K = filename[filename.find('K-')+2:filename.find('p-')]
    title = 'N=' + N + ' K=' + K

    data = np.loadtxt(filename, skiprows=2, delimiter=',')

    r = np.sqrt(data[:,0])
    n1 = data[:,1]
    n2 = data[:,2]
    t = data[:,3]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax.set_title(title)
    ax2.set_ylim(0,1)

    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Order parameter r')


    ax.plot(t, n1, 'r-')
    ax.plot(t, n2, 'g-')
    ax2.plot(t, r, 'b-')

plt.show()
