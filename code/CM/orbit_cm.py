import numpy as np
from pygadgetreader import *
import sys
import argparse

#parser = argparse.ArgumentParser(
#    description='Computing the CM position and velocity',
#    epilog='python orbit_cm.py snap_base_name initial_snap_number' \
#    'final_snap_number')

#output = parser.parse_args()


if len(sys.argv) != 4:
    print 'Usage: python orbit_cm.py snap_base_name inital_snap_number final_snap_number'
    print 'Ex: python orbit_cm.py snap 0 50'
    sys.exit(0)


# global variables to read the snapshots
snap = str(sys.argv[1])
# Initial and final snapshot number
i_n = float(sys.argv[2])
i_f = float(sys.argv[3])
path = '../../data/LMCMW/MW1LMC4/snapshots/'

# Number of Snapshots
N_snaps = (i_f - i_n) + 1
delta = 50 # precision of the CM computation in Kpc.
#Position and velocity arrays for the host and the satellite
X = np.zeros(N_snaps)
Y = np.zeros(N_snaps)
Z = np.zeros(N_snaps)

VX = np.zeros(N_snaps)
VY = np.zeros(N_snaps)
VZ = np.zeros(N_snaps)

Xsat = np.zeros(N_snaps)
Ysat = np.zeros(N_snaps)
Zsat = np.zeros(N_snaps)

VXsat = np.zeros(N_snaps)
VYsat = np.zeros(N_snaps)
VZsat = np.zeros(N_snaps)

#Galactocentric distance and time arrays
Rgal = np.zeros(N_snaps)
Vgal = np.zeros(N_snaps)
time = np.zeros(N_snaps)

# Defining function that computes the CM of the halo: 
def CM(x, y, z, delta):

    N = len(x) # Numero de particulas
    xCM = sum(x)/N
    yCM = sum(y)/N
    zCM = sum(z)/N

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
        # distance from the CM to all the particles
        R = np.sqrt((x - xCM)**2 + (y - yCM)**2 + (z - zCM)**2)
        # Finding the largest distance/velocity from the CM
        Rmax = max(R)
        # Selecting particles within half of the maximum radius
        index = np.where(r<Rmax/2.0)
        x = x[index]
        y = y[index]
        z = z[index]
        #Computing new CM
        xCM_new = sum(x)/len(x)
        yCM_new = sum(y)/len(y)
        zCM_new = sum(z)/len(z)
    return xCM_new, yCM_new, zCM_new

for i in range(0,len(X)):
    if i<10:
        time[i] = readheader(path + snap + "_00" + str(i),'time')
        positions = readsnap(path + snap + "_00" + str(i),'pos', 'dm')
        velocities = readsnap(path + snap + "_00" + str(i), 'vel', 'dm')
        particles_ids = readsnap(path + snap + "_00" + str(i), 'pid', 'dm')
    else:
        time[i] = readheader(path + snap + "_0" + str(i),'time')
        positions = readsnap(path + snap + "_0" + str(i),'pos', 'dm')
        velocities = readsnap(path + snap + "_0" + str(i), 'vel', 'dm')
        particles_ids = readsnap(path + snap + "_0" + str(i), 'pid', 'dm')

    X = np.sort(particles_ids)
    # The first half of particles are from the host DM halo, the 
    # second half are from the satellite DM halo.
    idcut = int(len(X)/2.0 - 1.0)
    index_mw = np.where(particles_ids<=X[idcut])
    index_LMC = np.where(particles_ids>X[idcut])

    x_mw = positions[index_mw[0],0]
    y_mw = positions[index_mw[0],1]
    z_mw = positions[index_mw[0],2]
    x_lmc = positions[index_LMC[0],0]
    y_lmc = positions[index_LMC[0],1]
    z_lmc = positions[index_LMC[0],2]

    vx_mw = velocities[index_mw[0],0]
    vy_mw = velocities[index_mw[0],1]
    vz_mw = velocities[index_mw[0],2]
    vx_lmc = velocities[index_LMC[0],0]
    vy_lmc = velocities[index_LMC[0],1]
    vz_lmc = velocities[index_LMC[0],2]

    X[i], Y[i], Z[i] = CM(x_mw, y_mw, z_mw, delta)
    Xsat[i], Ysat[i], Zsat[i]  = CM(x_lmc, y_lmc, z_lmc, delta)
    Rgal[i] = np.sqrt((X[i] - Xsat[i])**2 + (Y[i]-Ysat[i])**2 + (Z[i] - Zsat[i])**2)

    VX[i], VY[i], VZ[i] = CM(vx_mw, vy_mw, vz_mw, delta)
    VXsat[i], VYsat[i], VZsat[i]  = CM(vx_lmc, vy_lmc, vz_lmc, delta)
    Vgal[i] = np.sqrt((VX[i] - VXsat[i])**2 + (VY[i]-VYsat[i])**2 + (VZ[i] - VZsat[i])**2)


f = open("rgal_snaps_l0.txt", 'w')
f.write("#Time(Gyrs) | Rgal(kpc) | Xsat[kpc] | Ysat[kpc] | Zsat[kpc] |Xhost[kpc] | Yhost[kpc] Zhost[kpc] |" \
        " Vgal | Vxsat | Vysat | Vzsat | Vxhost | Vyhost | Vzhost | ")
for i in range(0, len(Rgal)):
    f.write("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"%(time[i], Rgal[i], Xsat[i], Ysat[i], \
    Zsat[i], X[i], Y[i], Z[i], VXsat[i], VYsat[i], VZsat[i], VX[i], VY[i], VZ[i]))
f.close()
