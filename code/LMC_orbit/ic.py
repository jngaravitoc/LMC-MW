import numpy as np 
from astropy import units
import sys

filename = sys.argv[1]
Rvir = float(sys.argv[2])
orbit = np.loadtxt(filename)

t = orbit[:,0]
x_LMC = orbit[:,1]
y_LMC = orbit[:,2]
z_LMC = orbit[:,3]
x_MW = orbit[:,4]
y_MW = orbit[:,5]
z_MW = orbit[:,6]
vx_LMC = orbit[:,7] 
vy_LMC = orbit[:,8]
vz_LMC = orbit[:,9]
vx_MW = orbit[:,10]
vy_MW = orbit[:,11]
vz_MW = orbit[:,12]


R = np.sqrt( (x_LMC - x_MW)**2 + (y_LMC - y_MW)**2 + (z_LMC - z_MW)**2 )
MIN = abs(R - Rvir)
index = np.where(MIN == min(MIN))

x_ic = x_LMC[index] - x_MW[index]
y_ic = y_LMC[index] - y_MW[index]
z_ic = z_LMC[index] - z_MW[index]
vx_ic = vx_LMC[index] - vx_MW[index]
vy_ic = vy_LMC[index] - vy_MW[index]
vz_ic = vz_LMC[index] - vz_MW[index]
vx_ic = vx_ic * units.kpc / units.Gyr
vy_ic = vy_ic * units.kpc / units.Gyr
vz_ic = vz_ic * units.kpc / units.Gyr

vx_ic = vx_ic.to(units.km / units.s)
vy_ic = vy_ic.to(units.km / units.s)
vz_ic = vz_ic.to(units.km / units.s)

print "x(kpc)", x_ic 
print "y(kpc)", y_ic
print "z(kpc)", z_ic
print "vx(km/s)", vx_ic
print "vy(km/s)", vy_ic
print "vz(km/s)", vz_ic


