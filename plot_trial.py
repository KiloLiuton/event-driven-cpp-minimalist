#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1];

data = np.loadtxt(filename, skiprows=2, delimiter=',')

r2 = data[:,0]
n1 = data[:,1]
n2 = data[:,2]

fig = plt.figure()
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax.plot(range(len(n1)), n1, 'r-')
ax.plot(range(len(n2)), n2, 'g-')
ax2.plot(range(len(r2)), r2, 'b-')

plt.show()
