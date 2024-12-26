#! /bin/python
import sys

import matplotlib

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    print(f"usage: {sys.argv[0]} file_in file_out")
    sys.exit()

y = np.loadtxt(sys.argv[1], dtype=float)

size = len(y)

x = np.linspace(0, 10, size)

#plt.ylim(-0.06, 0.06)
plt.ylim(1.0, 1.15)
plt.plot(x, y)
plt.savefig(sys.argv[2])

plt.show()
