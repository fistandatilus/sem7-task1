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

matplotlib.rcParams.update({'font.size': 15})
plt.ylim(-0.07, 0.07)
#plt.ylim(0.992,1.01)
plt.plot(x, y)
plt.savefig(sys.argv[2])

plt.show()
