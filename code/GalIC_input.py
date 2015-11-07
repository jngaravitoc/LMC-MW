import numpy as np 
import sys
from astropy import constants, units

Mvir = float(sys.argv[1])
Rvir = float(sys.argv[2])
Rs = float(sys.argv[3])


#Function that return the f function of the NFW profile.
def fx(x):
    f = np.log(1.+x) - (x / (1. + x))
    return f


# Function that computes the c200 for the NFW
def c(cvir, c200):
    q = 2.058
    y = (c200 / cvir) - (fx(c200) / (q * fx(cvir)))**(1./3.)
    return y


# Function to compute the c200
def bissection(cvir):
    min_c200 = 0.1
    max_c200 = cvir
    c_init = 0.5*(min_c200 + max_c200)
    y = c(cvir, c_init)
    while abs(y) > 0.000002:
        #print c_init
        if y>0:
            max_c200 = c_init
        if (y<0) :
            min_c200 = c_init
        c_init = 0.5*(min_c200 + max_c200)
        y = c(cvir, c_init)
    return c_init


def m200mvir(c200, cvir):
    x = fx(c200) / fx(cvir)
    return x

# function that computes the a/rs
def ars(c):
    x = 1 / ((2.0*fx(c))**(-0.5) - (1.0/c))
    return x

#Function that computes Mh/Mvir
def mhmvir(ar, cvir):
    x = ar**2 / (2.0*fx(cvir))
    return x


def v200(M):
    M = M * units.Msun
    G = constants.G
    H = 3.2407789E-18 * 0.7 / units.s
    v = (M * 10 * G * H)**(1.0/3.0)
    v = v.to(units.km/units.s)
    return v


def vvir(M):
    M = M * units.Msun
    G = constants.G
    H = 3.2407789E-18 * 0.7 / units.s
    v = (M*6.9713*G*H)**(1.0/3.0)
    v = v.to(units.km / units.s)
    return v


print "This code assume that your initial parameters are Mvir(NFW), Rvir(NFW) and Rs(NFW)"
cvir = Rvir / Rs
print "cvir(NFW) = ", cvir
c200 = bissection(cvir)
print "c200(NFW) = ", c200

a_rs_200 = ars(c200)
a_rs_vir = ars(cvir)
a200 = Rs * a_rs_vir
avir = Rs * a_rs_200

print "a200 = ", a200
print "avir = ", avir

M200 = Mvir * m200mvir(c200, cvir)

print "M200 = ", M200

mh_m200 = mhmvir(a_rs_200, c200)
mh_mvir = mhmvir(a_rs_vir, cvir)
mh_vir = Mvir * mh_mvir 
mh_200 = M200 * mh_m200

print "M_h_200 = ", mh_200
print "M_h_vir = ", mh_vir

v_h_200 = v200(mh_200)
v_h_vir = vvir(mh_vir)

print "v_h_200 = ", v_h_200
print "v_h_vir = ", v_h_vir
