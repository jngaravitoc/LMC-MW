import numpy as np 
from astropy import units
import astropy as apy
from profiles import *
from MW import *
import sys



# Implemented leapfrog algorithm to integrate orbits! 

def leapfrog(n_points, h, x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic):
      
	#n_points = 3000
	#h = 0.001
 	# Creating the arrays to collect the data in each step of the integration
        # the imput units should be in Kpc and Gyrs!

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

	#ax[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[0]
	#ay[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[1]
	#az[0] = acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[2]

        # one half step 

	t[1] = t[0] + h
	x[1] = x[0] + h * vx[0]
	y[1] = y[0] + h * vy[0]
	z[1] = z[0] + h * vz[0]

	vx[1] = vx[0] + h*acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[0]
	vy[1] = vy[0] + h*acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[1]
	vz[1] = vz[0] + h*acceleration(x[0], y[0], z[0], vx[0], vy[0], vz[0])[2]

	ax[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1])[0]
	ay[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1])[1]
	az[1] = acceleration(x[1],y[1], z[1], vx[1], vy[1], vz[1])[2]

	# iterate over all the steps!

	for i in range(2,n_points):
	    t[i] = t[i-1] + h
	    
	    x[i] = x[i-2] + 2 * h * vx[i-1]
	    y[i] = y[i-2] + 2 * h * vy[i-1]
	    z[i] = z[i-2] + 2 * h * vz[i-1]

	    vx[i] = vx[i-2] + 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1])[0]
	    vy[i] = vy[i-2] + 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1])[1]
	    vz[i] = vz[i-2] + 2 * h * acceleration(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1])[2]
	

        # Writing the output! 

	f = open('orbits.txt', 'w')
        f.write('#t(Gyrs), x(Kpc), y(Kpc), z(kpc), vx(kpc/gyr), vy(kpc/Gyr), vz(kpc/gyr) \n')
        for i in range(n_points):
		f.write(('%f \t %f \t %f \t %f \t %f \t %f \t %f \n')%(t[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]))
	f.close()

leapfrog(3000, 0.01, 40, 0, 0, 0, 215, 0)


