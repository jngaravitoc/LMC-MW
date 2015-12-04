import numpy as np 
from pygadgetreader import *
import sys

snap = str(sys.argv[1])
i_n = float(sys.argv[2])
i_f = float(sys.argv[3])

N_snaps = i_f - i_n + 1

X = np.zeros(N_snaps)
Y = np.zeros(N_snaps)
Z = np.zeros(N_snaps)

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

    while (np.sqrt((xCM_new-xCM)**2 + (yCM_new-yCM)**2
+(zCM_new-zCM)**2) > delta):
        xCM = xCM_new
        yCM = yCM_new
        zCM = zCM_new
        Rcm = np.sqrt(xCM**2 + yCM**2 + zCM**2)
        r = np.sqrt(x**2 + y**2 + z**2)
        R = np.sqrt((x - xCM)**2 + (y - yCM)**2 + (z - zCM)**2)
        Rmax = max(R)
        index = where(r<Rmax/2)
        x = x[index]
        y = y[index]
        z = z[index]
        #print Rmax
        xCM_new = sum(x)/len(x)
        yCM_new = sum(y)/len(y)
        zCM_new = sum(z)/len(z)
        #scatter(xCM_new, yCM_new)
    return xCM_new, yCM_new, zCM_new

for i in range(len(N_snaps)):
    positions = readsnap(snap + "_000",'pos', 'dm')
    particles_ids = readsnap("../../data/LMCMW/MW1LMC4/ICs1/MW1LMC4_000", 'pid', 'dm')
    x_sim_mw = positions[index_mw[0],0]
    y_sim_mw = positions[index_mw[0],1]
    z_sim_mw = positions[index_mw[0],2]
    x_sim_lmc = positions[index_LMC[0],0]
    y_sim_lmc = positions[index_LMC[0],1]
    z_sim_lmc = positions[index_LMC[0],2]
    X[i], Y[i], Z[i] = CM(x_sim_mw, y_sim_mw, z_sim_mw, 0, 0, 0, 10)

for i in range(len(X)):
    print X[i], Y[i], Z[i]
