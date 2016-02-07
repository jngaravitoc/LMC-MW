import numpy as np
from pygadgetreader import *
import sys

if len(sys.argv) != 5:
    print " python CM_resolution.py MW1LMC4a1H6_000 0.1 30000"
    sys.exit('Dude! put the parameters!')
snap = sys.argv[1]
deltar = float(sys.argv[2])
NhaloP = float(sys.argv[3]) # Number of particles in the halo.
nameout = sys.argv[4]

#Function that computes the CM of the halo, iteratively reducing the
#volume of the halo

def CM(x, y, z, vx, vy, vz, delta):
    N = len(x)
    xCM = sum(x)/len(x)
    yCM = sum(y)/len(y)
    zCM = sum(z)/len(z)

    xCM_new = np.zeros(100000)
    yCM_new = np.zeros(100000)
    zCM_new = np.zeros(100000)
    vxCM_new = np.zeros(100000)
    vyCM_new = np.zeros(100000)
    vzCM_new = np.zeros(100000)

    xCM_new[0] = xCM
    yCM_new[0] = yCM
    zCM_new[0] = zCM

    xCM = 0.0
    yCM = 0.0
    zCM = 0.0

    vxCM_new[0] = sum(vx)/N
    vyCM_new[0] = sum(vy)/N
    vzCM_new[0] = sum(vz)/N

    R1 = np.sqrt((x - xCM_new[0])**2 + (y - yCM_new[0])**2 + (z - zCM_new[0])**2)
    Rnow = [max(R1)]

    while (np.sqrt((xCM_new[-1]-xCM)**2 + (yCM_new[-1]-yCM)**2 +(zCM_new[-1]-zCM)**2) > delta):
        xCM = xCM_new[-1]
        yCM = yCM_new[-1]
        zCM = zCM_new[-1]
        Rcm = np.sqrt(xCM**2 + yCM**2 + zCM**2)
        R = np.sqrt((x - xCM)**2 + (y - yCM)**2 + (z - zCM)**2)
        Rmax = max(R)
        index = where(R<Rmax/2.0)
        x = x[index]
        y = y[index]
        z = z[index]
        vx = vx[index]
        vy = vy[index]
        vz = vz[index]
        N = len(x)
        xCM_new[i] = (sum(x)/N)
        yCM_new[i] = (sum(y)/N)
        zCM_new[i] = (sum(z)/N)
        vxCM_new[i] = (sum(vx)/N)
        vyCM_new[i] = (sum(vy)/N)
        vzCM_new[i] = (sum(vz)/N)
        Rnow.append(max(np.sqrt((x - xCM_new[-1])**2 + (y - yCM_new[-1])**2 + (z - zCM_new[-1])**2)))

    return xCM_new, yCM_new, zCM_new, vxCM_new, vyCM_new, vxCM_new, Rnow

#function that reads the snapshot

def loading_data(filename):
    pothalos = readsnap(filename, 'pot', 'dm')
    poshalos = readsnap(filename, 'pos', 'dm')
    velhalos = readsnap(filename, 'vel', 'dm')
    idhalos = readsnap(filename, 'pid', 'dm')


    potdisk = readsnap(filename, 'pot', 'disk')
    posdisk = readsnap(filename, 'pos', 'disk')
    veldisk = readsnap(filename, 'vel', 'disk')

    return pothalos, poshalos[:,0], poshalos[:,1], poshalos[:,2], \
           velhalos[:,0],  velhalos[:,1],  velhalos[:,2], idhalos, potdisk,\
           posdisk[:,0], posdisk[:,1], posdisk[:,2],  veldisk[:,0], veldisk[:,1],\
           veldisk[:,2]
# function that returns the MW and the LMC DM particles

def LMCMWparticles(ids, NMW, x, y, z, vx, vy, vz, pot):
    X = np.sort(ids)
    limit = X[NMW]
    index_MW = np.where(ids<=limit)[0]
    index_LMC = np.where(ids>limit)[0]
    xmw, ymw, zmw = x[index_MW], y[index_MW], z[index_MW]
    vxmw, vymw, vzmw = vx[index_MW], vy[index_MW], vz[index_MW]
    xlmc, ylmc, zlmc = x[index_LMC], y[index_LMC], z[index_LMC]
    vxlmc, vylmc, vzlmc = vx[index_LMC], vy[index_LMC], vz[index_LMC]
    potmw = pot[index_MW]
    potlmc = pot[index_LMC]
    return xmw, ymw, zmw, vxmw, vymw, vzmw, xlmc, ylmc, zlmc, vxlmc, vylmc, vzlmc, potmw, potlmc

# Function that returns the coordiantes of the particle with the minimum potential

def potentialCM(potential, x, y, z):
    minpotential = np.where(potential == min(potential))
    XCM, YCM, ZCM = x[minpotential], y[minpotential], z[minpotential]
    return XCM, YCM, ZCM

Ph, Xh, Yh, Zh, VxH, VyH, VzH, idH, Pd, Xd, Yd, Zd, Vxd, Vyd, Vzd = loading_data(snap)
XMW, YMW, ZMW, VxMW, VyMW, VzMW, XL, YL, ZL, VxL, VyL, VzL, PMW, PL = LMCMWparticles(idH, NhaloP, Xh, Yh, Zh, VxH, VyH, VzH, Ph)

# Computing CM with the potential
XCMPMW, YCMPMW, ZCMPMW = potentialCM(PMW, XMW, YMW, ZMW)
XCMPd, YCMPd, ZCMPd = potentialCM(Pd, Xd, Yd, Zd)
XCMPL, YCMPL, ZCMPL = potentialCM(PL, XL, YL, ZL)

# Shell model
XCMMW, YCMMW, ZCMMW, vXCMMW, vYCMMW, vZCMMW, RsMW = CM(XMW, YMW, ZMW, VxMW, VyMW, VzMW, deltar)
XCMD, YCMD, ZCMD, vXCMD, vYCMD, vZCMD, RsD = CM(Xd, Yd, Zd, Vxd, Vyd, Vzd, deltar)
XCML, YCML, ZCML, vXCML, vYCML, vZCML, RsL = CM(XL, YL, ZL, VxL, VyL, VzL, deltar)

# Computing galactocentric radius
Rmw = np.sqrt(XCMMW**2 + YCMMW**2 + ZCMMW**2)
Vmw = np.sqrt(vXCMMW**2 + vYCMMW**2 + vZCMMW**2)

Rlmc = np.sqrt(XCMMW**2 + YCMMW**2 + ZCMMW**2)
Vlmc = np.sqrt(vXCMMW**2 + vYCMMW**2 + vZCMMW**2)

Rd = np.sqrt(XCMD**2 + YCMD**2 + ZCMD**2)
Vd = np.sqrt(vXCMD**2 + vYCMD**2 + vZCMD**2)

RCMPMW = np.sqrt(XCMPMW**2 + YCMPMW**2 + ZCMPMW**2)
RCMPD = np.sqrt(XCMPd**2 + YCMPd**2 + ZCMPd**2)
RCMPL = np.sqrt(XCMPL**2 + YCMPL**2 + ZCMPL**2)

fmw = open('MW'+nameout, 'w')
fmw.write('RShell(kpc), Rmw, Vmw, Rpot \n')
for i in range(len(RsMW)):
    fmw.write(('%f %f %f %f\n')%(RsMW[i], Rmw[i], Vmw[i], RCMPMW))
fmw.close()

fd = open('Disk' + nameout, 'w')
fd.write('RShell(kpc), Rdisk, Vdisk, RDpot \n')
for i in range(len(RsD)):
    fd.write(('%f %f %f %f\n')%(RsD[i], Rd[i], Vd[i], RCMPMW))
fmw.close()


fl = open('LMC' + nameout, 'w')
fl.write('RShell(kpc), Rdisk, Vdisk, RDpot \n')
for i in range(len(RsL)):
    fl.write(('%f %f %f %f\n')%(RsL[i], Rlmc[i], Vlmc[i], RCMPL))
fmw.close()

