import numpy as np 
import matplotlib.pyplot as plt
import sys


filename = sys.argv[1]

data = np.loadtxt(filename)

x = data[:,0]
y = data[:,1]
z = data[:,2]
R = data[:,3]
t = data[:,4]



plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.plot(x, z)
plt.xlabel('x')
plt.ylabel('z')
plt.show()

plt.plot(t, R)
plt.xlabel('t')
plt.ylabel('R')
plt.show()
