import numpy as np
from pygadgetreader import *
import sys

if len(sys.argv) != 5:
    print "Ex:"
    print "python CM_resolution.py MW1LMC4a1H6_000 0.1 30000 outname"
    sys.exit('Error: Not enough parameters')
snap = sys.argv[1]
deltar = float(sys.argv[2])
NhaloP = float(sys.argv[3]) # Number of particles in the halo.
nameout = sys.argv[4]

#Function that computes the CM of the halo, iteratively reducing the
#volume of the halo

def CM(x, y, z, vx, vy, vz, delta):
    N = len(x)
    xCM = sum(x)/N
    yCM = sum(y)/N
    zCM = sum(z)/N

    xCM_new = np.zeros(300)
    yCM_new = np.zeros(300)
    zCM_new = np.zeros(300)
    vxCM_new = np.zeros(300)
    vyCM_new = np.zeros(300)
    vzCM_new = np.zeros(300)
    Rnow = np.zeros(300)

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
    Rnow[0] = max(R1)
    i=0
    while (np.sqrt((xCM_new[i]-xCM)**2 + (yCM_new[i]-yCM)**2 +(zCM_new[i]-zCM)**2) > delta):
        xCM = xCM_new[i]
        yCM = yCM_new[i]
        zCM = zCM_new[i]
        Rcm = np.sqrt(xCM**2 + yCM**2 + zCM**2)
        R = np.sqrt((x - xCM)**2 + (y - yCM)**2 + (z - zCM)**2)
        Rmax = max(R)
        index = np.where(R<Rmax/1.3)[0]
        x = x[index]
        y = y[index]
        z = z[index]
        vx = vx[index]
        vy = vy[index]
        vz = vz[index]
        N = len(x)
        i+=1
        xCM_new[i] = (sum(x)/N)
        yCM_new[i] = (sum(y)/N)
        zCM_new[i] = (sum(z)/N)
        vxCM_new[i] = (sum(vx)/N)
        vyCM_new[i] = (sum(vy)/N)
        vzCM_new[i] = (sum(vz)/N)
        Rnow[i] = max(np.sqrt((x - xCM_new[i])**2 + (y - yCM_new[i])**2 + (z - zCM_new[i])**2))
    clean = np.where(Rnow != 0)[0]
    return xCM_new[clean], yCM_new[clean], zCM_new[clean], vxCM_new[clean], vyCM_new[clean], vzCM_new[clean], Rnow[clean]

def potential_CM(potential, x, y, z, vx, vy, vz):
    index = np.where(potential< min(potential)*0.90)[0]
    x_p = x[index]
    y_p = y[index]
    z_p = z[index]
    vx_p = vx[index]
    vy_p = vy[index]
    vz_p = vz[index]
    N = len(x_p)
    print N
    x_cm = sum(x_p)/N
    y_cm = sum(y_p)/N
    z_cm = sum(z_p)/N
    vx_cm = sum(vx_p)/N
    vy_cm = sum(vy_p)/N
    vz_cm = sum(vz_p)/N
    Rcm = np.sqrt(x_cm**2.0 + y_cm**2.0 + z_cm**2.0)
    Vcm = np.sqrt(vx_cm**2.0 + vy_cm**2.0 + vz_cm**2.0)
    return x_cm, y_cm, z_cm, vx_cm, vy_cm, vz_cm, Rcm, Vcm

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
    index_MW = np.where(ids<limit)[0]
    index_LMC = np.where(ids>=limit)[0]
    xmw, ymw, zmw = x[index_MW], y[index_MW], z[index_MW]
    vxmw, vymw, vzmw = vx[index_MW], vy[index_MW], vz[index_MW]
    xlmc, ylmc, zlmc = x[index_LMC], y[index_LMC], z[index_LMC]
    vxlmc, vylmc, vzlmc = vx[index_LMC], vy[index_LMC], vz[index_LMC]
    potmw = pot[index_MW]
    potlmc = pot[index_LMC]
    return xmw, ymw, zmw, vxmw, vymw, vzmw, xlmc, ylmc, zlmc, vxlmc, vylmc, vzlmc, potmw, potlmc

# Function that returns the coordiantes of the particle with the minimum potential

def potentialCM(potential, x, y, z):
    #rvir = np.where(np.sqrt(x**2 + y**2 + z**2) < 200)[0]
    #potential, x, y, z = potential[rvir], x[rvir], y[rvir], z[rvir]
    minpotential = np.where(potential == min(potential))
    XCM, YCM, ZCM = x[minpotential], y[minpotential], z[minpotential]
    return XCM, YCM, ZCM

def potentialLCM(potential, x, y, z, Xcm, Ycm, Zcm):
    r = np.sqrt((x-Xcm)**2 + (y-Ycm)**2 + (z-Zcm)**2)
    rvir = np.where(r<100)
    potential, x, y, z = potential[rvir], x[rvir], y[rvir], z[rvir]
    minpotential = np.where(potential == min(potential))
    XCM, YCM, ZCM = x[minpotential], y[minpotential], z[minpotential]
    return XCM, YCM, ZCM

Ph, Xh, Yh, Zh, VxH, VyH, VzH, idH, Pd, Xd, Yd, Zd, Vxd, Vyd, Vzd = loading_data(snap)
print 'Number of DM particles: ', len(Xh)
print 'NUmber of Disk particles: ', len(Xd)
XMW, YMW, ZMW, VxMW, VyMW, VzMW, XL, YL, ZL, VxL, VyL, VzL, PMW, PL = LMCMWparticles(idH, NhaloP, Xh, Yh, Zh, VxH, VyH, VzH, Ph)

# Shell model
XCMMW, YCMMW, ZCMMW, vXCMMW, vYCMMW, vZCMMW, RsMW = CM(XMW, YMW, ZMW, VxMW, VyMW, VzMW, deltar)
XCMD, YCMD, ZCMD, vXCMD, vYCMD, vZCMD, RsD = CM(Xd, Yd, Zd, Vxd, Vyd, Vzd, deltar)
XCML, YCML, ZCML, vXCML, vYCML, vZCML, RsL = CM(XL, YL, ZL, VxL, VyL, VzL, deltar)

# Computing CM with the potential
XCMPMW, YCMPMW, ZCMPMW = potentialCM(PMW, XMW, YMW, ZMW, VxMW, VyMW, VzMW)
#XCMPd, YCMPd, ZCMPd = potentialCM(Pd, Xd, Yd, Zd)
XCMPL, YCMPL, ZCMPL = potentialCM(PL, XL, YL, ZL, XCML[-1], YCML[-1], ZCML[-1])

# Computing galactocentric radius
Rmw = np.sqrt(XCMMW**2 + YCMMW**2 + ZCMMW**2)
Vmw = np.sqrt(vXCMMW**2 + vYCMMW**2 + vZCMMW**2)

Rlmc = np.sqrt(XCML**2 + YCML**2 + ZCML**2)
Vlmc = np.sqrt(vXCML**2 + vYCML**2 + vZCML**2)

Rd = np.sqrt(XCMD**2 + YCMD**2 + ZCMD**2)
Vd = np.sqrt(vXCMD**2 + vYCMD**2 + vZCMD**2)

RCMPMW = np.sqrt(XCMPMW**2 + YCMPMW**2 + ZCMPMW**2)
RCMPD = np.sqrt(XCMPd**2 + YCMPd**2 + ZCMPd**2)
RCMPL = np.sqrt(XCMPL**2 + YCMPL**2 + ZCMPL**2)

#------------------- Writing data ---------------------

fmw = open('MW'+nameout, 'w')
fmw.write('RShell(kpc), Rmw, Vmw, Rpot \n')
for i in range(len(RsMW)):
    fmw.write(('%f %f %f %f \n')%(RsMW[i], Rmw[i], Vmw[i], RCMPMW))
fmw.close()

fd = open('Disk'+nameout, 'w')
fd.write('RShell(kpc), Rdisk, Vdisk, RDpot \n')
for i in range(len(RsD)):
    fd.write(('%f %f %f %f \n')%(RsD[i], Rd[i], Vd[i], RCMPD))
fd.close()

print len(RsL), len(Rlmc), len(Vlmc), len(RCMPL)
fl = open('LMC'+nameout, 'w')
fl.write('RShell(kpc), Rlmc, Vlmc, RLMCpot \n')
for i in range(len(RsL)):
    fl.write(('%f %f %f %f \n')%(RsL[i], Rlmc[i], Vlmc[i], RCMPL))
fl.close()
