import numpy as np 
from profiles import *
from cosmotools import *
from math import erf
import sys

x_ic  =  -1
y_ic  =  -41
z_ic  =  -28
vx_ic  =  -57
vy_ic  =  -226
vz_ic  =  221
#x_ic = 10
#y_ic = 0
#z_ic = 0
#vx_ic = -10
#vy_ic = 0
#vz_ic = 0
Mbulge = 1E10
Mhalo = float(sys.argv[1])
Mdisk = float(sys.argv[2])
Msat = float(sys.argv[3])
Rvir = float(sys.argv[4])  
rs = float(sys.argv[5]) * units.kpc
rlmc = float(sys.argv[6]) * units.kpc
mw = sys.argv[7]

ra = 3.5 
rb = 0.53

print '#Initial conditions:'
print '#x = ', x_ic, '(kpc),',  'y = ', y_ic, '(kpc),',  'z = ', z_ic, '(kpc),',  'vx = ', vx_ic, '(km/s),' , 'vy = ', vy_ic,  '(km/s),', 'vz = ', vz_ic, '(km/s),'


def coulomb_log(r):
    bmax = r # position of test particle at a time t
    k = 3 * units.kpc # kpc
    bmin = 1.4 * k # k is the softening length if the LMC were modeled using a plummer progile . See Besla07
    L = bmax / bmin
    return np.log(L)

def sigma(rs, r, M_halo, Rv):
    M_halo = M_halo * units.Msun
    Rvir = Rv * units.kpc 
    vvir = np.sqrt( G * M_halo / Rvir) 
    c = Rvir.value / rs.value
    g = np.log(1+c) - (c /(1+c))
    vmax = np.sqrt(0.216 * vvir**2 * c / g)
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
    Mhalo = M_halo - M_disk - M_bulge
    c = Rvir / rs.value
    rho = dens_NFWnRvir(c, x.value, y.value, z.value, Mhalo, Rvir) #+ dens_mn(ra, rb, x.value, y.value, z.value, M_disk) + dens_hernquist(0.7, r.value, M_bulge) # a, r, v  ****
    # Mass of the satellite
    M_sat = M_sat * units.Msun
    # Computing the dyanmical friction
    G = constants.G
    G = G.to(units.kpc**3 / units.Msun / units.Gyr**2)
    factor = - 4 * np.pi * G**2
    Coulomb =  coulomb_log(r) #**************
    s = sigma(rs, r, M_halo, Rvir) 
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
    c = Rvir / rs.value
    ahalo = a_NFWnRvir(c, x, y, z, M_halo, Rvir)
    adisk = a_mn(ra, rb, x, y, z, M_disk)
    abulge = a_hernquist(0.7, x, y, z, M_bulge)
    ax = ahalo[0] + adisk[0] + abulge[0]
    ay = ahalo[1] + adisk[1] + abulge[1]
    az = ahalo[2] + adisk[2] + abulge[2]
    ax = ax.to(units.kpc/units.Gyr**2)  
    ay = ay.to(units.kpc/units.Gyr**2) 
    az = az.to(units.kpc/units.Gyr**2)
    r = np.sqrt(x**2 + y**2 + z**2)

    # Truncate the halo at the virial radius

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


def acceleration_mw(x, y, z, M_sat):
    x = x * units.kpc
    y = y * units.kpc
    z = z * units.kpc
    M_sat = M_sat * units.Msun
    r = np.sqrt(x**2 + y**2 + z**2)
    G = constants.G
    #rlmc = rlmc * units.kpc
    #if (r.value<=Rvir):
    Ax =  - G * M_sat * x / (r**2 + rlmc**2)**(3.0/2.0)        
    Ay =  - G * M_sat * y / (r**2 + rlmc**2)**(3.0/2.0)
    Az =  - G * M_sat * z / (r**2 + rlmc**2)**(3.0/2.0)
    #else: 
#	Ax = - G * M_sat * x / r**3
#       Ay = - G * M_sat * y / r**3
#       Az = - G * M_sat * z / r**3
	
    Ax = Ax.to(units.kpc / units.Gyr**2)
    Ay = Ay.to(units.kpc / units.Gyr**2)
    Az = Az.to(units.kpc / units.Gyr**2)
    return Ax.value, Ay.value, Az.value	


def leapfrog(x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic, M_halo, M_disk, M_bulge, M_sat, Rvir):

        n_points = 8000
        h = 0.001

        # Creating the arrays to collect the data in each step of the integration
        # the imput units should be in Kpc for positions and km/s for velocities
 
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
        x = np.zeros(n_points) #+ x_ic
        y = np.zeros(n_points) #+ y_ic
        z = np.zeros(n_points) #+ z_ic

        x_mw = np.zeros(n_points)
        y_mw = np.zeros(n_points)
        z_mw = np.zeros(n_points)


        vx = np.zeros(n_points)
        vy = np.zeros(n_points)
        vz = np.zeros(n_points)

        vx_mw  = np.zeros(n_points)
        vy_mw  = np.zeros(n_points)
        vz_mw  = np.zeros(n_points)


        ax = np.zeros(n_points)
        ay = np.zeros(n_points)
        az = np.zeros(n_points)

        ax_mw = np.zeros(n_points)
        ay_mw = np.zeros(n_points)
        az_mw = np.zeros(n_points)

        t[0] = 0

        # This initial conditions come form MW.py, the units are Kpc and Gyr
        x[0] = x_ic
        y[0] = y_ic
        z[0] = z_ic
        
        x_mw[0] = 0.0
        y_mw[0] = 0.0
        z_mw[0] = 0.0

        vx[0] = vx_ic
        vy[0] = vy_ic
        vz[0] = vz_ic

        vx_mw[0] = 0.0
        vy_mw[0] = 0.0
        vz_mw[0] = 0.0 

        ax[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
        ay[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
        az[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

       	ax_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), M_sat)[0]
        ay_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), M_sat)[1]
        az_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), M_sat)[2]


        # one half step 
        # at the first step the MW is at 0, 0, 0 and its not moving
 
        t[1] = t[0] - h
        x[1] = x[0] - h * vx[0]
        y[1] = y[0] - h * vy[0]
        z[1] = z[0] - h * vz[0]
        
        if (mw=='free'):
        	x_mw[1] = x_mw[0] - h * vx_mw[0]
        	y_mw[1] = y_mw[0] - h * vy_mw[0]
        	z_mw[1] = z_mw[0] - h * vz_mw[0]

		vx_mw[1] = vx_mw[0] - h * ax_mw[0]
        	vy_mw[1] = vy_mw[0] - h * ay_mw[0]
        	vz_mw[1] = vz_mw[0] - h * az_mw[0]

		ax_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[0]
        	ay_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[1]
        	az_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[2]


        vx[1] = vx[0] - h * ax[0]
        vy[1] = vy[0] - h * ay[0]
        vz[1] = vz[0] - h * az[0]

        ax[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
        ay[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
        az[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

        #ax_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[0] 
        #ay_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[1]
        #az_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), M_sat)[2]

        # iterate over all the steps!

        for i in range(2,n_points):
            t[i] = t[i-1] - h

            x[i] = x[i-2] - 2 * h * vx[i-1] 
            y[i] = y[i-2] - 2 * h * vy[i-1] 
            z[i] = z[i-2] - 2 * h * vz[i-1] 

            if (mw=='free'):
            	x_mw[i] = x_mw[i-2] - 2 * h * vx_mw[i-1]
            	y_mw[i] = y_mw[i-2] - 2 * h * vy_mw[i-1]
            	z_mw[i] = z_mw[i-2] - 2 * h * vz_mw[i-1]
        
                vx_mw[i] = vx_mw[i-2] - 2 * h * ax_mw[i-1]
                vy_mw[i] = vy_mw[i-2] - 2 * h * ay_mw[i-1]
                vz_mw[i] = vz_mw[i-2] - 2 * h * az_mw[i-1]

	        ax_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[0]
                ay_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[1]
                az_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[2]


            vx[i] = vx[i-2] - 2 * h * ax[i-1]
            vy[i] = vy[i-2] - 2 * h * ay[i-1]
            vz[i] = vz[i-2] - 2 * h * az[i-1]
      
            
            #vx_mw[i] = vx_mw[i-2] - 2 * h * ax_mw[i-1]
            #vy_mw[i] = vy_mw[i-2] - 2 * h * ay_mw[i-1]
            #vz_mw[i] = vz_mw[i-2] - 2 * h * az_mw[i-1]

	    ax[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
	   
            ay[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
	    
            az[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]
	
	    #ax_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[0]
            #ay_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[1]
            #az_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), M_sat)[2]

        return t, x, y, z, x_mw, y_mw, z_mw, vx, vy, vz, vx_mw, vy_mw, vz_mw, ax, ay, az, ax_mw, ay_mw, az_mw


t, xx , yy, zz, xx_mw, yy_mw, zz_mw, vx, vy, vz, vx_mw, vy_mw, vz_mw, ax, ay, az, ax_mw, ay_mw, az_mw = leapfrog(x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic, Mhalo, Mbulge, Mdisk, Msat, Rvir)

#R = np.sqrt((xx-xx_mw)**2 + (yy-yy_mw)**2 + (zz-zz_mw)**2)
R = np.sqrt((xx)**2 + (yy)**2 + (zz)**2)
RMW = np.sqrt(xx_mw**2 + yy_mw**2 + zz_mw**2)

print "# t, x, y, z, xmw, ymw, zmw, vxmw, vymw, vzmw, ax, ay, az, axmw, aymw, azmw"
for i in range(len(xx)):
	print t[i], xx[i], yy[i], zz[i], xx_mw[i], yy_mw[i], zz_mw[i], vx[i], vy[i], vz[i],  vx_mw[i], vy_mw[i], vz_mw[i], ax[i], ay[i], az[i], ax_mw[i], ay_mw[i], az_mw[i] 


#print "# t, x, y, z, R"
#for i in range(len(xx)):
#	print t[i], xx_mw[i], yy_mw[i], zz_mw[i], RMW[i] 

