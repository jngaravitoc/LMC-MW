"""
A short python code to compute the SCF coefficients.

Nicolas Garavito-Camargo, University of Arizona.

08/10/16:
Implementation of the phi_nl_f, Anl_f and the coefficients functions.
note: The code reproduce the results from biff. github.com/adrn/biff

"""

import numpy as np
from scipy import special
from math import factorial


def phi_nl_f(r, n, l):
    """
    Input:
    ------
    r: Particles coordinates
    n: n
    l: l
    Output:
    -------
    phi_nl: The potential basis. Eq.11 Lowing e.a (2011)
    """

    factor = r**l * (1.+r)**(-2.*l-1.) * np.sqrt(4.*np.pi) 
    s = (r-1.)/(r+1.)
    C_n = special.eval_gegenbauer(n, 2.*l+3./2., s)
    return -factor*C_n.real

def Anl_f(n,l):
    """
    Computes Eq 16. Lowing e.a (2011)
    """
    K_nl = 0.5*n*(n + 4.*l + 3.) + (l + 1.)*(2.*l + 1.)
    factor = 2.**(8.*l + 6.)/(4.*np.pi*K_nl)
    return -factor*factorial(n)*(n + 2*l+ 3./2.)*(special.gamma(2.*l +
3./2.))**2.0/special.gamma(n + 4.*l + 3.)


def coefficients(x, y, z, mass, n, l, m):
    """
    Input:
    ------
    x: Array of particles x-coordiante
    y: Array of particles y-coordiante
    z: Array of particles z-coordinate
    mass: Array with the mass particles
    n: n
    l: l
    m: m
    r_s: halo scale length
    Output:
    -------
    var(a): A float number with the variance
    of the coefficient (S_nlm, T_nlm) evaluated at n,l,m.
    """

    r = np.sqrt(x**2.0+y**2.0+z**2.0)/r_s
    N = len(r)
    theta = np.arccos(z/(r*r_s))
    phi = np.arctan2(y,x)
    Y_lm = special.sph_harm(m,l,0,theta) # scipy notation m,l
    phi_nl = phi_nl_f(r, n, l)
    Psi = phi_nl*Y_lm*mass
    if m==0:
        dm0=1.0
    else:
        dm0=0.0
    Anl = Anl_f(n,l)
    Snlm=(2.-dm0)*Anl*np.sum(Psi.real*np.cos(m*phi))
    Tnlm=(2.-dm0)*Anl*np.sum(Psi.real*np.sin(m*phi))

    return Snlm, Tnlm
