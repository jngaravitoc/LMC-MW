import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


filename = sys.argv[1]

data = np.loadtxt(filename)

t = data[:,0]
x = data[:,1]
y = data[:,2]
z = data[:,3]
xmw = data[:,4]
ymw = data[:,5]
zmw = data[:,6]


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z)
plt.scatter(x[0], y[0], z[0], c='r')
ax.plot(xmw, ymw, zmw)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


Rsat = np.sqrt(x**2 + y**2 + z**2)
Rmw = np.sqrt(xmw**2 + ymw**2 + zmw**2)
Rgal = np.sqrt((x-xmw)**2 + (y-ymw)**2 + (z-zmw)**2)
plt.plot(t, Rsat)
plt.plot(t, Rmw)
plt.plot(t, Rgal)
plt.xlabel('t')
plt.ylabel('R')
plt.ylim(0, 1000)
plt.show()
