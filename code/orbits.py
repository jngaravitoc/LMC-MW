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
ax.scatter(x[-1], y[-1], z[-1])
ax.plot(xmw, ymw, zmw)
ax.scatter(xmw[-1], ymw[-1], zmw[-1])
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('Sagorbit.png', bbox_inches='tight')
plt.show()


plt.plot(x-xmw, z-zmw)
plt.xlabel('$x$', fontsize=25)
plt.ylabel('$z$', fontsize=25)
plt.scatter(x[-1], z[-1])
plt.savefig('SagXZ.png', bbox_inches='tight')
plt.show()


Rsat = np.sqrt(x**2 + y**2 + z**2)
Rmw = np.sqrt(xmw**2 + ymw**2 + zmw**2)
Rgal = np.sqrt((x-xmw)**2 + (y-ymw)**2 + (z-zmw)**2)
#plt.plot(t, Rsat)
#plt.plot(t, Rmw)
plt.plot(t, Rgal)
plt.xlabel('t(Gyr)')
plt.ylabel('Rgal(Kpc)')
#plt.axhline(261)
#plt.ylim(0, 1000)
plt.savefig('Saggalradius.png', bbox_inches='tight')
plt.show()
