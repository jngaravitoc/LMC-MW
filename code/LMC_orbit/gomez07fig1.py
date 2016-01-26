import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3]
filename4 = sys.argv[4]
filename5 = sys.argv[5]
filename6 = sys.argv[6]


data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)
data5 = np.loadtxt(filename5)
data6 = np.loadtxt(filename6)

t = data1[:,0]
x1= data1[:,1]
y1= data1[:,2]
z1= data1[:,3]
xmw1 = data1[:,4]
ymw1 = data1[:,5]
zmw1 = data1[:,6]

x2= data2[:,1]
y2= data2[:,2]
z2= data2[:,3]
xmw2 = data2[:,4]
ymw2 = data2[:,5]
zmw2 = data2[:,6]


x3= data3[:,1]
y3= data3[:,2]
z3= data3[:,3]
xmw3 = data3[:,4]
ymw3 = data3[:,5]
zmw3 = data3[:,6]

x4= data4[:,1]
y4= data4[:,2]
z4= data4[:,3]
xmw4 = data4[:,4]
ymw4 = data4[:,5]
zmw4 = data4[:,6]

x5= data5[:,1]
y5= data5[:,2]
z5= data5[:,3]
xmw5 = data5[:,4]
ymw5 = data5[:,5]
zmw5 = data5[:,6]

x6= data6[:,1]
y6= data6[:,2]
z6= data6[:,3]
xmw6 = data6[:,4]
ymw6 = data6[:,5]
zmw6 = data6[:,6]



Rsat1 = np.sqrt(x1**2 + y1**2 + z1**2)
Rmw1 = np.sqrt(xmw1**2 + ymw1**2 + zmw1**2)
Rgal1 = np.sqrt((x1-xmw1)**2 + (y1-ymw1)**2 + (z1-zmw1)**2)

Rsat2 = np.sqrt(x2**2 + y2**2 + z2**2)
Rmw2 = np.sqrt(xmw2**2 + ymw2**2 + zmw2**2)
Rgal2 = np.sqrt((x2-xmw2)**2 + (y2-ymw2)**2 + (z2-zmw2)**2)

Rsat3 = np.sqrt(x3**2 + y3**2 + z3**2)
Rmw3 = np.sqrt(xmw3**2 + ymw3**2 + zmw3**2)
Rgal3 = np.sqrt((x3-xmw3)**2 + (y3-ymw3)**2 + (z3-zmw3)**2)

Rsat4 = np.sqrt(x4**2 + y4**2 + z4**2)
Rmw4 = np.sqrt(xmw4**2 + ymw4**2 + zmw4**2)
Rgal4 = np.sqrt((x4-xmw4)**2 + (y4-ymw4)**2 + (z4-zmw4)**2)

Rsat5 = np.sqrt(x5**2 + y5**2 + z5**2)
Rmw5 = np.sqrt(xmw5**2 + ymw5**2 + zmw5**2)
Rgal5 = np.sqrt((x5-xmw5)**2 + (y5-ymw5)**2 + (z5-zmw5)**2)

Rsat6 = np.sqrt(x6**2 + y6**2 + z6**2)
Rmw6 = np.sqrt(xmw6**2 + ymw6**2 + zmw6**2)
Rgal6 = np.sqrt((x6-xmw6)**2 + (y6-ymw6)**2 + (z6-zmw6)**2)
#plt.plot(t, Rsat)
#plt.plot(t, Rmw)

plt.plot(t, Rgal1, ls='-', c='b', lw=2)
plt.plot(t, Rgal2, ls='-', c='r', lw=2)
plt.plot(t, Rgal3, ls='-', c='g', lw=2)
plt.plot(t, Rgal4, ls='-', c='k', lw=2)
plt.plot(t, Rgal5, ls='-', c='m', lw=2)
plt.plot(t, Rgal6, ls='-', c='y', lw=2)
plt.axhline(329, c='k', alpha=0.7)
plt.xlabel('$\mathrm{Look\ Back\ Time(Gyr)}$', fontsize=22)
plt.ylabel('$\mathrm{Galactocentric\ Radius(kpc)}$', fontsize=22)
plt.ylim(min(Rgal1), 750)
plt.savefig("gomez.png", bbox_inches='tight')
plt.show()
