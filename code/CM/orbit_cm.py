import numpy as np
from pygadgetreader import *
import sys
import argparse

if len(sys.argv) != 8:
    print 'Usage: python orbit_cm.py snap_base_name inital_snap_number'\
           'final_snap_number path2file out?name  #DMhost #DMsat'
    print 'Ex: python orbit_cm.py snap 0 50 pat2file out_name Nhost Nsat'
    sys.exit(0)


# global variables to read the snapshots
snap = str(sys.argv[1])
# Initial and final snapshot number
i_n = int(sys.argv[2])
i_f = int(sys.argv[3])
out_name = str(sys.argv[5])
Nhost = int(sys.argv[6])
Nsat = int(sys.argv[7])
path = str(sys.argv[4])#'../../data/LMCMW/MW1LMC4/a1/'

# Number of Snapshots
N_snaps = (i_f - i_n) + 1

deltar = 0.5 # precision of the CM computation in Kpc.

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

#Function that computes the CM of the halo:
def CMMW(x, y, z, pot):
    rcut = np.where(np.sqrt(x**2+y**2+z**2)<50)[0]
    x, y, z, pot = x[rcut], y[rcut], z[rcut], pot[rcut]
    cm = np.where(pot == min(pot))[0]
    x_cm, y_cm, z_cm = x[cm], y[cm], z[cm]
    return x_cm, y_cm, z_cm

def CMLMC(x, y, z, pot, xcmmw, ycmmw, zcmmw):
    rcut = np.where(np.sqrt((x-xcmmw)**2+(y-ycmmw)**2+(z-zcmmw)**2)<10)[0]
    x, y, z, pot = x[rcut], y[rcut], z[rcut], pot[rcut]
    cm = np.where(pot == min(pot))[0]
    x_cm, y_cm, z_cm = x[cm], y[cm], z[cm]
    return x_cm, y_cm, z_cm

def VCM(x, y, z, xcm, ycm, zcm, vx, vy, vz):
    Ntot = len(x)
    N = Ntot
    while(N>0.1*Ntot):
        rshell = np.sqrt((x-xcm)**2 + (y-ycm)**2 + (z-zcm)**2)
        rcut = rshell / 1.2
        cut = np.where(rshell<=rcut)[0]
        x, y, z = x[cut], y[cut], z[cut]
        vx, vy, vz = vx[cut], vy[cut], vz[cut]
        N = len(X)
    vxcm = sum(vx)/N
    vycm = sum(vy)/N
    vzcm = sum(vz)/N
    return vxcm, vycm, vzcm


for i in range(i_n, i_f + 1):
    if i<10:
        time[i-i_n] = readheader(path + snap + "_00" + str(i),'time')
        positions = readsnap(path + snap + "_00" + str(i),'pos', 'dm')
        velocities = readsnap(path + snap + "_00" + str(i), 'vel', 'dm')
        particles_ids = readsnap(path + snap + "_00" + str(i), 'pid', 'dm')
        potential = readsnap(path + snap + "_00" + str(i), 'pid', 'dm')
    elif ((i>=10) & (i<100)):
        time[i-i_n] = readheader(path + snap + "_0" + str(i),'time')
        positions = readsnap(path + snap + "_0" + str(i),'pos', 'dm')
        velocities = readsnap(path + snap + "_0" + str(i), 'vel', 'dm')
        particles_ids = readsnap(path + snap + "_0" + str(i), 'pid', 'dm')
        potential = readsnap(path + snap + "_0" + str(i), 'pid', 'dm')
    else:
        time[i-i_n] = readheader(path + snap + "_" + str(i),'time')
        positions = readsnap(path + snap + "_" + str(i),'pos', 'dm')
        velocities = readsnap(path + snap + "_" + str(i), 'vel', 'dm')
        particles_ids = readsnap(path + snap + "_" + str(i), 'pid', 'dm')
        potential = readsnap(path + snap + "_" + str(i), 'pid', 'dm')

    ID = np.sort(particles_ids)
    # The first set of particles are from the host DM halo, the
    # second set are from the satellite DM halo, the limit is know by
    # the number of particles in the host halo.
    idcut = ID[Nhost-1]
    index_mw = np.where(particles_ids<=idcut)
    index_LMC = np.where(particles_ids>idcut)

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

    potmw = potential[index_mw]
    potlmc = potential[index_LMC]

    X[i-i_n], Y[i-i_n], Z[i-i_n] =  CMMW(x_mw, y_mw, z_mw, potmw)
    Xsat[i-i_n], Ysat[i-i_n], Zsat[i-i_n] = CMLMC(x_lmc, y_lmc, z_lmc, potlmc, X[i-i_n], Y[i-i_n], Z[i-i_n])
    VX[i-i_n], VY[i-i_n], VZ[i-i_n] = VCM(x_mw, y_mw, z_mw, X, Y, Z, vx_mw, vy_mw, vz_mw)
    VXsat[i-i_n], VYsat[i-i_n], VZsat[i-i_n] = VCM(x_lmc, y_lmc, z_lmc, Xsat, Ysat, Zsat, vx_lmc, vy_lmc, vz_lmc)

    Rgal[i-i_n] = np.sqrt((X[i-i_n] - Xsat[i-i_n])**2 + (Y[i-i_n]-Ysat[i-i_n])**2 + (Z[i-i_n] - Zsat[i-i_n])**2)
    Vgal[i-i_n] = np.sqrt((VX[i-i_n] - VXsat[i-i_n])**2 + (VY[i-i_n]-VYsat[i-i_n])**2 + (VZ[i-i_n] - VZsat[i-i_n])**2)
    print Rgal, Vgal, X, Y, Z, Xsat, Ysat, Zsat, VX, VY, VZ, VXsat, VYsat, VZsat
"""
f = open(out_name, 'w')
f.write("#Time(Gyrs) | Rgal(kpc) | Xsat[kpc] | Ysat[kpc] | Zsat[kpc] |Xhost[kpc] | Yhost[kpc] Zhost[kpc] |"\
        "Vgal | Vxsat | Vysat | Vzsat | Vxhost | Vyhost | Vzhost |\n")

for i in range(0, len(Rgal)):
    f.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(time[i], Rgal[i], Xsat[i], Ysat[i],\
    Zsat[i], X[i], Y[i], Z[i], Vgal[i], VXsat[i], VYsat[i], VZsat[i], VX[i], VY[i], VZ[i]))
f.close()
"""
