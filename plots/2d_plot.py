#! /bin/python
import sys

import matplotlib

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 4:
    print(f"usage: {sys.argv[0]} T_st file_in file_out")
    sys.exit()
z = np.loadtxt(sys.argv[2], dtype=float)

x_size = len(z[0])
y_size = len(z)

x = np.linspace(0, 10, x_size+1)
y = np.linspace(0, float(sys.argv[1]), y_size+1)


#Z = z.reshape(y_size, x_size - 1)  

X, Y = np.meshgrid(x, y)


plt.pcolormesh(Y, X, z, shading='flat', cmap='rainbow')
plt.ylabel('Пространство')
plt.xlabel('Время')
plt.colorbar()
plt.savefig(sys.argv[3])

plt.show()
