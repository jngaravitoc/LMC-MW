import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pygadgetreader import *
import sys

#function to read data

if len(sys.argv)!=6:
    sys.exit('Error: wrong number of parameters')
path = sys.argv[1]
snap = sys.argv[2]
NparticlesMW = float(sys.argv[3])
Ni = float(sys.argv[4])
Nf = float(sys.argv[5])
N = int(Nf - Ni + 1.0)

def loading_data(filename):
    pothalos = readsnap(filename, 'pot', 'dm')
    poshalos = readsnap(filename, 'pos', 'dm')
    velhalos = readsnap(filename, 'vel', 'dm')
    idhalos = readsnap(filename, 'pid', 'dm')


    potdisk = readsnap(filename, 'pot', 'disk')
    posdisk = readsnap(filename, 'pos', 'disk')
    veldisk = readsnap(filename, 'vel', 'disk')

    return poshalos[:,0], poshalos[:,1], poshalos[:,2], idhalos
# function that returns the MW and the LMC DM particles

def LMCMWparticles(ids, NMW, x, y, z):
    X = np.sort(ids)
    limit = X[NMW]
    index_MW = np.where(ids<limit)[0]
    index_LMC = np.where(ids>=limit)[0]
    xmw, ymw, zmw = x[index_MW], y[index_MW], z[index_MW]
    xlmc, ylmc, zlmc = x[index_LMC], y[index_LMC], z[index_LMC]
    return xmw, ymw, zmw, xlmc, ylmc, zlmc

def plot(xmw, ymw, zmw, xlmc, ylmc, zlmc, snap):
    index = np.where(np.sqrt(xmw**2 + ymw**2 + zmw**2)<300.0)[0]
    plt.figure(figsize=(16, 5))
    plt.subplot(1, 3, 1)
    plt.hist2d(xmw[index], ymw[index], bins=200, norm=mpl.colors.LogNorm(), cmap=plt.get_cmap('hot'))
    plt.subplot(1, 3, 2)
    plt.hist2d(xmw[index], zmw[index], bins=200, norm=mpl.colors.LogNorm(), cmap=plt.get_cmap('hot'))
    plt.subplot(1, 3, 3)
    plt.hist2d(ymw[index], zmw[index], bins=200, norm=mpl.colors.LogNorm(), cmap=plt.get_cmap('hot'))
    plt.savefig(snap + '.png', bbox_inches='tight')

for i in range(N):
    if i<10:
        xDM, yDM, zDM, idDM = loading_data(path+snap+ '_00'+str(i))
    if i>+10:
        xDM, yDM, zDM, idDM = loading_data(path+snap+ '_0'+str(i))

    xMW, yMW, zMW, xLMC, yLMC, zLMC = LMCMWparticles(idDM, NparticlesMW, xDM, yDM, zDM)
    plot(xMW, yMW, zMW, xLMC, yLMC, zLMC, snap+str(i))
