import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


filename1 = sys.argv[1]
filename2 = sys.argv[2]

data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)


t1 = data1[:,0]
x1 = data1[:,1]
y1 = data1[:,2]
z1 = data1[:,3]
xmw1 = data1[:,4]
ymw1 = data1[:,5]
zmw1 = data1[:,6]

t2 = data2[:,0]
x2 = data2[:,1]
y2 = data2[:,2]
z2 = data2[:,3]
xmw2 = data2[:,4]
ymw2 = data2[:,5]
zmw2 = data2[:,6]



#Rsat = np.sqrt(x**2 + y**2 + z**2)
#Rmw = np.sqrt(xmw**2 + ymw**2 + zmw**2)
Rgal1 = np.sqrt((x1-xmw1)**2 + (y1-ymw1)**2 + (z1-zmw1)**2)
Rgal2 = np.sqrt((x2-xmw2)**2 + (y2-ymw2)**2 + (z2-zmw2)**2)
#plt.plot(t, Rsat)
#plt.plot(t, Rmw)
plt.plot(t1, Rgal1, ls='--', c='b', lw=2 )
plt.plot(t2, Rgal2, ls='-', c='b', lw=2)
plt.xlabel('t(Gyr)')
plt.ylabel('Rgal(Kpc)')
plt.axhline(261)
#plt.ylim(0, 1000)
plt.savefig('Sag_orbits.png', bbox_inches='tight')
plt.show()
