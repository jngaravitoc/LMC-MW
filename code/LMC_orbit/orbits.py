import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


filename1 = sys.argv[1]
data1 = np.loadtxt(filename1)

filename2 = sys.argv[2]
data2 = np.loadtxt(filename2)

filename3 = sys.argv[3]
data3 = np.loadtxt(filename3)

filename4 = sys.argv[4]
data4 = np.loadtxt(filename4)

filename5 = sys.argv[5]
data5 = np.loadtxt(filename5)

filename6 = sys.argv[6]
data6 = np.loadtxt(filename6)


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

t3 = data3[:,0]
x3 = data3[:,1]
y3 = data3[:,2]
z3 = data3[:,3]
xmw3 = data3[:,4]
ymw3 = data3[:,5]
zmw3 = data3[:,6]

t4 = data4[:,0]
x4 = data4[:,1]
y4 = data4[:,2]
z4 = data4[:,3]
xmw4 = data4[:,4]
ymw4 = data4[:,5]
zmw4 = data4[:,6]

t5 = data5[:,0]
x5 = data5[:,1]
y5 = data5[:,2]
z5 = data5[:,3]
xmw5 = data5[:,4]
ymw5 = data5[:,5]
zmw5 = data5[:,6]

t6 = data6[:,0]
x6 = data6[:,1]
y6 = data6[:,2]
z6 = data6[:,3]
xmw6 = data6[:,4]
ymw6 = data6[:,5]
zmw6 = data6[:,6]


#Rsat = np.sqrt(x**2 + y**2 + z**2)
#Rmw = np.sqrt(xmw**2 + ymw**2 + zmw**2)
Rgal1 = np.sqrt((x1-xmw1)**2 + (y1-ymw1)**2 + (z1-zmw1)**2)
Rgal2 = np.sqrt((x2-xmw2)**2 + (y2-ymw2)**2 + (z2-zmw2)**2)
Rgal3 = np.sqrt((x3-xmw3)**2 + (y3-ymw3)**2 + (z3-zmw3)**2)
Rgal4 = np.sqrt((x4-xmw4)**2 + (y4-ymw4)**2 + (z4-zmw4)**2)
Rgal5 = np.sqrt((x5-xmw5)**2 + (y5-ymw5)**2 + (z5-zmw5)**2)
Rgal6 = np.sqrt((x6-xmw6)**2 + (y6-ymw6)**2 + (z6-zmw6)**2)


#plt.plot(t, Rsat)
#plt.plot(t, Rmw)
plt.plot(t1, Rgal1, label='model1', lw=2)
plt.plot(t2, Rgal2, label='model2', lw=2)
plt.plot(t3, Rgal3, label='model3', lw=2)
plt.plot(t4, Rgal4, label='model4', lw=2)
plt.plot(t5, Rgal5, label='model5', lw=2)
plt.plot(t6, Rgal6, label='model6', lw=2)
plt.xlabel('t(Gyr)')
plt.ylabel('Rgal(Kpc)')
plt.legend(loc='best',fontsize=15)
plt.axhline(261)
#plt.ylim(0, 1000)
plt.savefig('LMC_orbits.png', bbox_inches='tight')
plt.show()
