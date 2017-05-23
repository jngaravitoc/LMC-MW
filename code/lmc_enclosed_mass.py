 # -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pygadgetreader import *
import sys

def MW_LMC_particles(xyz, mass, pids, NMW_particles):
    """
    Function that return the MW and the LMC particles
    positions and velocities.

    Parameters:
    -----------
    xyz: snapshot coordinates with shape (n,3)
    vxys: snapshot velocities with shape (n,3)
    pids: particles ids
    NMW_particles: Number of MW particles in the snapshot
    Returns:
    --------
    xyz_mw, vxyz_mw, xyzlmc, vxyz_lmc: coordinates and velocities of
    the MW and the LMC.

    """
    sort_indexes = np.sort(pids)
    N_cut = sort_indexes[NMW_particles]
    MW_ids = np.where(pids<N_cut)[0]
    LMC_ids = np.where(pids>=N_cut)[0]
    return xyz[LMC_ids], mass[LMC_ids]


def re_center(pos, cm):
    """
    Re center a halo to its center of mass.
    """
    pos_n = np.copy(pos)
    for i in range(3):
        pos_n[:,i] = pos[:,i] - cm[i]
    return pos_n


def enclosed_mass(pos, mass):
    """
    Computes the enclosed mass of a halo
    
    """
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    dr = np.linspace(0.5, 15, 20)
    cm = np.zeros(len(dr))
    for i in range(len(dr)):
        index = np.where(r<dr[i])
        cm[i] = np.sum(mass[index])
    return dr, cm

if __name__ == "__main__":

    path = sys.argv[1]
    snap_base = sys.argv[2]
    orbit = np.loadtxt('LMC_orbit/orbits/LMC6_40Mb0_orbit.txt')
    # r = np.zeros((30,20))
    # cm_r = np.zeros((30, 20))
    for i in range(1, 150, 30):
        mwlmc_pos = readsnap(path+snap_base+'{:0>3d}'.format(i), 'pos', 'dm')
        mwlmc_mass = readsnap(path+snap_base+'{:0>3d}'.format(i), 'mass', 'dm')
        mwlmc_pid = readsnap(path+snap_base+'{:0>3d}'.format(i), 'pid', 'dm')
        time = readhead(path+snap_base+'{:0>3d}'.format(i), 'time')
        orbit = np.loadtxt('LMC_orbit/orbits/LMC6_40Mb0_orbit.txt')
        lmc_pos, lmc_mass = MW_LMC_particles(mwlmc_pos, mwlmc_mass, mwlmc_pid, 37500000)
        lmc_pos_cm = re_center(lmc_pos, orbit[i][6:9])
        r, cm_r = enclosed_mass(lmc_pos_cm, lmc_mass)
        plt.semilogy(r, cm_r*1E10, label='{:.1f}'.format(time))
    plt.axhline(1.7E10, ls='--', c='k', alpha=0.5)
    plt.axvline(8.7, ls='--', c='k', alpha=0.5)
    plt.legend()
    plt.show()
