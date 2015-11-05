import numpy as np 
import sys
from astropy import constants
from astropy import units

print 'input: Mvir(Msun), Cvir :'
Mvir = float(sys.argv[1])
cvir = float(sys.argv[2])
print Mvir,'Msun', cvir
def fx(x):
    f = np.log(1.+x) - ( x / (1. + x) )
    return f

def ars(c):
    x = 1 / ( (2.0*fx(c))**(-0.5) - (1.0/c) )
    return x

def mhmvir(ar, cvir):
    x = ar**2 / (2.0*fx(cvir))
    return x




a_rs = ars(cvir)

mh_mvir = mhmvir(a_rs, cvir)

mh = mh_mvir * Mvir
mh = mh * units.Msun
G = constants.G
H = 3.2407789E-18 * 0.7 / units.s
v = (mh * 10 * G * H)**(1.0/3.0) 
v = v.to(units.km/units.s)
v_h = (Mvir*units.Msun *10*G*H)**(1.0/3.0)
v_h = v_h.to(units.km / units.s)

Rvir = v/(H*10)
A = Rvir / cvir
A = A.to(units.kpc)


print 'Hernquist velocity:'
print v_h


print '-------------In case you want the Hernquist equivalence of a NFW Mvir -----------'

print 'This is the velocity you need = ', v
print 'This is the mass of Hernquist equivalent to a NFW', mh
print 'Hernquist scale length = ', A
