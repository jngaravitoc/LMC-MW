"""
Implementation of the algorithm to compute the coefficients
of the SCF method. Following the Appendix Weinberg 1996.

History:
--------

06/17/2016:
Implemented 
"""

import numpy as np
from scipy import linalg

def smooth_b(a):
    b = np.zeros(len(a))
    for i in range(len(a)):
        b[i] = 1/(1.0 + np.var(a)/a[i]**2)
    return b

def MISE(a):
    b = smooth_b(a)
    D = np.sum(b**2 * np.var(a) + (b-1)**2.0*a**2.0)
    return D

def S_matrix(a):
    a = a.flatten()
    S = np.zeros((len(a), len(a)))
    for i in range(len(a)):
        for j in range(len(a)):
            S[i][j] = a[i]*a[j]
    return S

def lambda_prime(lambdas):
    F = np.zeros(len(lambdas))
    print len(F)
    for i in range(len(F)):
        F[i] = np.sum(lambdas.real[:i])/np.sum(lambdas.real)
    return F


def a_prime(T, a):
    a_new = np.zeros(len(a))
    for i in range(len(a_new)):
        a_new[i] = np.sum(np.dot(T.real[i],a))
    return a_new

