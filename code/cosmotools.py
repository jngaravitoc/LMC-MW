import numpy as np
from astropy import units
from astropy import constants

H_0 = 70 * units.km / (units.s * units.Mpc)
Omega0 = 0.27
G = constants.G
G = G.to(units.kiloparsec**3 / (units.Msun * units.s**2))

def H(z):
    Lambda0 = 1. - Omega0
    return H_0*(Omega0*(1+z)**3 - (Omega0+Lambda0-1)*(1+z)**2 + Lambda0)**0.5

def Omega_z(z):
    return Omega0 * (1+z)**3 * (H_0/H(z))**2

def rho_crit(z):
    H2 = H(z)**2
    rho = 3*H2 / (8*np.pi*G)
    return rho

def Dvir(z):# from the solution of the top hat model! 
    Omegaz = Omega_z(z)
    x = Omegaz - 1
    Deltavir =  ( (18*np.pi**2) +  (82*x) - 39*x**2 ) / Omegaz
    return Deltavir

def rvir(Mvir, z):
    Mvir = Mvir * units.Msun
    Deltavir = Dvir(z)
    pcrit = rho_crit(z)
    Rvir = ( 3*Mvir / (4 * np.pi * Deltavir * pcrit * Omega0) )**(1/3.)
    Rvir = Rvir.to(units.kpc)
    return Rvir



