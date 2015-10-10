import numpy as np 
from astropy import units
import astropy as apy
from profiles import *
from math import erf

G = apy.constants.G
G = G.to(units.kpc**3 / units.Msun / units.s**2)


def coulomb_log(r):
    bmax = r # position of test particle at a time t
    k = 3 * units.kpc # kpc
    bmin = 1.6 * k # k is the softening length if the LMC were modeled using a plummer progile . See Besla07
    L = bmax / bmin
    return np.log(L)

def dynamical_friction_sis(x, y, z, vx, vy, vz, M_sat):
    x = x * units.kpc
    y = y * units.kpc
    z = z * units.kpc
    r = np.sqrt(x**2 + y**2 + z**2)
    vx = vx * units.km / units.s
    vy = vy * units.km / units.s
    vz = vz * units.km / units.s
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    a = 10 # concentration parameter
    # Density at distance r and a velocity v
    rho = dens_sis(10, r.value, v.value) # a, r, v
    v = v.to(units.kpc / units.s)
    M_sat = M_sat * units.Msun
    factor = - 4 * np.pi * G**2  
    Coulomb = coulomb_log(r)
    sigma = v / np.sqrt(2)
    X = v / ( np.sqrt(2) * sigma ) #Check this
    F_dfx = factor * M_sat * rho * Coulomb / v**3 * (  erf(X) - 2*X/(np.sqrt(np.pi) * np.exp(-X**2))  ) * vx
    F_dfy = factor * M_sat * rho * Coulomb / v**3 * (  erf(X) - 2*X/(np.sqrt(np.pi) * np.exp(-X**2))  ) * vy
    F_dfz = factor * M_sat * rho * Coulomb / v**3 * (  erf(X) - 2*X/(np.sqrt(np.pi) * np.exp(-X**2))  ) * vz
    F_dfx = F_dfx.to(units.kpc / units.Gyr**2)
    F_dfy = F_dfy.to(units.kpc / units.Gyr**2)
    F_dfz = F_dfz.to(units.kpc / units.Gyr**2)
    return F_dfx.value, F_dfy.value, F_dfz.value


def acceleration(x, y, z, vx, vy, vz):
    M_bulge = 1E10
    M_disk = 5.5E10
    M_halo = 1E12
    M_sat = 1E11
    abulge = a_hernquist(0.7, x, y, z, M_bulge)
    adisk = a_mn(6.5, 0.6, x, y, z, M_disk)
    ahalo = a_NFW(11.0, x, y, z, M_halo)
    #print abulge, adisk, ahalo
    ax = abulge[0] + adisk[0] + ahalo[0]
    ay = abulge[1] + adisk[1] + ahalo[1]
    az = abulge[2] + adisk[2] + ahalo[2] 
    ax = ax.to(units.kpc/units.Gyr**2) # / G.value / M_bulge
    ay = ay.to(units.kpc/units.Gyr**2) #/ G.value / M_bulge
    az = az.to(units.kpc/units.Gyr**2) #/ G.value / M_bulge
    F_dfx, F_dfy, F_dfz = dynamical_friction_sis(x, y, z, vx, vy, vz, M_sat)
    Ax = ax.value - F_dfx
    Ay = ay.value - F_dfy
    Az = az.value - F_dfz
    return Ax, Ay, Az

