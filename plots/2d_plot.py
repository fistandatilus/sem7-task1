#! /bin/python
import sys

import matplotlib

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 5:
    print(f"usage: {sys.argv[0]} x_mesh_size y_mesh_size file T_st")
    sys.exit()

x_size = int(sys.argv[1])
y_size = int(sys.argv[2])

x = np.linspace(0, 10, x_size+1)
y = np.linspace(0, float(sys.argv[4]), y_size+1)

z = np.loadtxt(sys.argv[3], dtype=float)

Z = z.reshape(y_size, x_size - 1)  

X, Y = np.meshgrid(x[:-1], y)


plt.pcolormesh(Y, X, Z, shading='flat', cmap='rainbow')
plt.ylabel('Пространство')
plt.xlabel('Время')
plt.colorbar()
plt.savefig('z_V_.jpg')

plt.show()
