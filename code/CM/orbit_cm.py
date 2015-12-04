import numpy as np 
from pygadgetreader import *
import sys

snap = str(sys.argv[1])
i_n = float(sys.argv[2])
i_f = float(sys.argv[3])
path = '/home/xozidok/work/github/LMC-MW/data/LMCMW/MW1LMC4/snapshots/'


N_snaps = (i_f - i_n) + 1

X = np.zeros(N_snaps)
Y = np.zeros(N_snaps)
Z = np.zeros(N_snaps)

Xsat = np.zeros(N_snaps)
Ysat = np.zeros(N_snaps)
Zsat = np.zeros(N_snaps)

Rgal = np.zeros(N_snaps)
time = np.zeros(N_snaps)

def CM(x, y, z, xCM, yCM, zCM, delta):
    xCM = sum(x)/len(x)
    yCM = sum(y)/len(y)
    zCM = sum(z)/len(z)

    xCM_new = xCM
    yCM_new = yCM
    zCM_new = zCM

    xCM = 0.0
    yCM = 0.0
    zCM = 0.0

    while ((np.sqrt((xCM_new-xCM)**2 + (yCM_new-yCM)**2 \
          +(zCM_new-zCM)**2) > delta)):
        xCM = xCM_new
        yCM = yCM_new
        zCM = zCM_new
        Rcm = np.sqrt(xCM**2 + yCM**2 + zCM**2)
        r = np.sqrt(x**2 + y**2 + z**2)
        R = np.sqrt((x - xCM)**2 + (y - yCM)**2 + (z - zCM)**2)
        Rmax = max(R)
        index = np.where(r<Rmax/2)
        x = x[index]
        y = y[index]
        z = z[index]
        xCM_new = sum(x)/len(x)
        yCM_new = sum(y)/len(y)
        zCM_new = sum(z)/len(z)
    return xCM_new, yCM_new, zCM_new

for i in range(0,len(X)):
    if i<10:
        time[i] = readheader(path + snap + "_00" + str(i),'time')
        positions = readsnap(path + snap + "_00" + str(i),'pos', 'dm')
        particles_ids = readsnap(path + snap + "_00" + str(i), 'pid', 'dm')
    else:
        time[i] = readheader(path + snap + "_0" + str(i),'time')
        positions = readsnap(path + snap + "_0" + str(i),'pos', 'dm')
        particles_ids = readsnap(path + snap + "_0" + str(i), 'pid', 'dm')
    X = np.sort(particles_ids)
    index_mw = np.where(particles_ids<=49376)
    index_LMC = np.where(particles_ids>49376)
    x_mw = positions[index_mw[0],0]
    y_mw = positions[index_mw[0],1]
    z_mw = positions[index_mw[0],2]
    x_lmc = positions[index_LMC[0],0]
    y_lmc = positions[index_LMC[0],1]
    z_lmc = positions[index_LMC[0],2]
    X[i], Y[i], Z[i] = CM(x_mw, y_mw, z_mw, 0, 0, 0, 100)
    Xsat[i], Ysat[i], Zsat[i]  = CM(x_lmc, y_lmc, z_lmc, 0, 0, 0, 100)
    Rgal[i] = np.sqrt((X[i] - Xsat[i])**2 + (Y[i]-Ysat[i])**2 + (Z[i] - Zsat[i])**2)

f = open("rgal_snaps_l0.txt", 'w')
for i in range(0, len(Rgal)):
    f.write("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"%(time[i], Rgal[i], Xsat[i], Ysat[i], \
    Zsat[i], X[i], Y[i], Z[i]))
f.close()
