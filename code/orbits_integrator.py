import numpy as np 
from profiles import *
from cosmotools import *
from math import erf
import sys

x_ic  =  float(sys.argv[1])
y_ic  =  float(sys.argv[2])
z_ic  =  float(sys.argv[3])
vx_ic  =  float(sys.argv[4])
vy_ic  =  float(sys.argv[5])
vz_ic  =  float(sys.argv[6])
Mhalo = float(sys.argv[7])
Mdisk = float(sys.argv[8])
Mbulge = float(sys.argv[9])
Msat = float(sys.argv[10])
Rvir = 261

print '#Initial conditions:'
print '#x = ', x_ic, '(kpc),',  'y = ', y_ic, '(kpc),',  'z = ', z_ic, '(kpc),',  'vx = ', vx_ic, '(km/s),' , 'vy = ', vy_ic,  '(km/s),', 'vz = ', vz_ic, '(km/s),'


def coulomb_log(r):
    bmax = r # position of test particle at a time t
    k = 3 * units.kpc # kpc
    bmin = 1.4 * k # k is the softening length if the LMC were modeled using a plummer progile . See Besla07
    L = bmax / bmin
    return np.log(L)

def sigma(c, r, M_halo, Rv):
    M_halo = M_halo * units.Msun
    Rvir = Rv * units.kpc #rvir2(M_halo.value, 0)  # Halo mass, Z
    vvir = np.sqrt( G * M_halo / Rvir) 
    g = np.log(1+c) - (c /(1+c))
    vmax = np.sqrt(0.216 * vvir**2 * c / g)
    rs = Rvir / c
    x = r / rs
    sigma = vmax * 1.4393 * x **(0.354) / (1 + 1.1756*x**0.725)
    sigma = sigma.to(units.kpc / units.Gyr)
    return sigma

def dynamical_friction_sis(x, y, z, vx, vy, vz, M_halo, M_disk, M_bulge, M_sat, Rvir):
    # Coordinates
    x = x * units.kpc
    y = y * units.kpc
    z = z * units.kpc
    r = np.sqrt(x**2 + y**2 + z**2)
    # Velocities
    vx = vx * units.kpc / units.Gyr
    vy = vy * units.kpc / units.Gyr
    vz = vz * units.kpc / units.Gyr
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    # Density of the NFW at a given r
    Mhalo = (M_halo - M_disk - M_bulge)
    rho = dens_NFWnRvir(11, x.value, y.value, z.value, M_halo, Rvir) + dens_mn(3.5, 0.53, x.value, y.value, z.value, M_disk) + dens_hernquist(0.7, r.value, M_bulge) # a, r, v  ****
    # Mass of the satellite
    M_sat = M_sat * units.Msun
    # Computing the dyanmical friction
    G = constants.G
    G = G.to(units.kpc**3 / units.Msun / units.Gyr**2)
    factor = - 4 * np.pi * G**2
    Coulomb =  coulomb_log(r) #**************
    s = sigma(11, r, M_halo, Rvir) 
    X = v / ( np.sqrt(2) * s ) 
    #print v, s, X
    # Main equation
    a_dfx = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vx) / v**3.0
    a_dfy = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vy) / v**3.0
    a_dfz = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vz) / v**3.0 
    # Transforming to the right units
    a_dfx = a_dfx.to(units.kpc / units.Gyr**2) 
    a_dfy = a_dfy.to(units.kpc / units.Gyr**2)
    a_dfz = a_dfz.to(units.kpc / units.Gyr**2)
    return a_dfx.value, a_dfy.value, a_dfz.value

def acceleration(x, y, z, vx, vy, vz, M_halo, M_disk, M_bulge, M_sat, Rvir):
    ahalo = a_NFWnRvir(11.0, x, y, z, M_halo, Rvir)
    adisk = a_mn(3.5, 0.53, x, y, z, M_disk)
    abulge = a_hernquist(0.7, x, y, z, M_bulge)
    ax = ahalo[0] + adisk[0] + abulge[0]
    ay = ahalo[1] + adisk[1] + abulge[1]
    az = ahalo[2] + adisk[2] + abulge[2]
    ax = ax.to(units.kpc/units.Gyr**2)  
    ay = ay.to(units.kpc/units.Gyr**2) 
    az = az.to(units.kpc/units.Gyr**2)
    r = np.sqrt(x**2 + y**2 + z**2)
    if (r<=Rvir):
    	a_dfx, a_dfy, a_dfz = dynamical_friction_sis(x, y, z, vx, vy, vz, M_halo, M_disk, M_bulge, M_sat, Rvir)
    	Ax = ax.value + a_dfx 
    	Ay = ay.value + a_dfy 
    	Az = az.value + a_dfz
    else:
        G = constants.G
        Mtot = (M_halo) * units.Msun
    	Ax = - G * Mtot * x * units.kpc / (r*units.kpc)**3 
        Ay = - G * Mtot * y * units.kpc / (r*units.kpc)**3 
        Az = - G * Mtot * z * units.kpc / (r*units.kpc)**3
        Ax = Ax.to(units.kpc / units.Gyr**2)
        Ay = Ay.to(units.kpc / units.Gyr**2)
        Az = Az.to(units.kpc / units.Gyr**2)
        Ax = Ax.value
        Ay = Ay.value
        Az = Az.value
    return Ax, Ay, Az  

def leapfrog(x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic, M_halo, M_disk, M_bulge, M_sat, Rvir):

        n_points = 8000
        h = 0.001

        # Creating the arrays to collect the data in each step of the integration
        # the imput units should be in Kpc and Gyrs!
 
        # Coordinates tranformation

        vx_ic = vx_ic * units.km / units.s
        vy_ic = vy_ic * units.km / units.s
        vz_ic = vz_ic * units.km / units.s
        vx_ic = vx_ic.to(units.kpc / units.Gyr) 
        vy_ic = vy_ic.to(units.kpc / units.Gyr)
        vz_ic = vz_ic.to(units.kpc / units.Gyr)
        vx_ic = vx_ic.value
        vy_ic = vy_ic.value
        vz_ic = vz_ic.value
        
        t = np.zeros(n_points)
        x = np.zeros(n_points)
        y = np.zeros(n_points)
        z = np.zeros(n_points)

        vx = np.zeros(n_points)
        vy = np.zeros(n_points)
        vz = np.zeros(n_points)


        ax = np.zeros(n_points)
        ay = np.zeros(n_points)
        az = np.zeros(n_points)

        t[0] = 0

        # This initial conditions come form MW.py, the units are Kpc and Gyr
        x[0] = x_ic
        y[0] = y_ic
        z[0] = z_ic


        vx[0] = vx_ic
        vy[0] = vy_ic
        vz[0] = vz_ic

        ax[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
        ay[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
        az[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

        # one half step 

        t[1] = t[0] - h
        x[1] = x[0] - h * vx[0]
        y[1] = y[0] - h * vy[0]
        z[1] = z[0] - h * vz[0]

        vx[1] = vx[0] - h * ax[0]
        vy[1] = vy[0] - h * ay[0]
        vz[1] = vz[0] - h * az[0]

        ax[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
        ay[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
        az[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

        # iterate over all the steps!

        for i in range(2,n_points):
            t[i] = t[i-1] - h

            x[i] = x[i-2] - 2 * h * vx[i-1]
            y[i] = y[i-2] - 2 * h * vy[i-1]
            z[i] = z[i-2] - 2 * h * vz[i-1]

            vx[i] = vx[i-2] - 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
            vy[i] = vy[i-2] - 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
            vz[i] = vz[i-2] - 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]
      
        
        return x, y, z, t,  


xx , yy, zz, t = leapfrog(x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic, Mhalo, Mbulge, Mdisk, Msat, Rvir)

R = np.sqrt(xx**2 + yy**2 + zz**2)

for i in range(len(xx)):
	print xx[i], yy[i], zz[i], R[i], t[i]
 
