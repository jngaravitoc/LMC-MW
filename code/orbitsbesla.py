import numpy as np 
import matplotlib.pyplot as plt
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

#x = data[:,0]
#y = data[:,1]
#z = data[:,2]
R1 = data1[:,3]
R2 = data2[:,3]
R3 = data3[:,3]
R4 = data4[:,3]
R5 = data5[:,3]
R6 = data6[:,3]

t = data1[:,4]


"""
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.plot(x, z)
plt.xlabel('x')
plt.ylabel('z')
plt.show()
"""
plt.plot(t, R1, ls='--', c='b', lw=2)
plt.plot(t, R2, ls='--', c='#FF8000', lw=2)
plt.plot(t, R3, ls='--', c='g', lw=2)
plt.plot(t, R4, ls='--', c='k', lw=2)
plt.plot(t, R5, ls='--', c='m', lw=2)
plt.plot(t, R6, ls='--', c='y', lw=2)
plt.axhline(261, c='k', alpha=0.7)
plt.xlabel('$\mathrm{Look\ Back\ Time(Gyr)}$', fontsize=22)
plt.ylabel('$\mathrm{Galactocentric\ Radius(kpc)}$', fontsize=22)
plt.ylim(min(R1), 1000)
plt.savefig('dfLMC.png', bbox_inches='tight')
plt.show()
