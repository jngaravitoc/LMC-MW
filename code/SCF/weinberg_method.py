"""
Implementation of the algorithm to compute the coefficients
of the SCF method. Following the Appendix Weinberg 1996.

History:
--------

06/17/2016:
Implemented the MISE method.
07/11/2016
Implementing the PCA method.
"""

import numpy as np
from scipy import linalg

# ----------------------- Optimal convergence MISE -------------------

def smooth_b(a):
    # Eq 8, Weinberg 96
    b = np.zeros(len(a))
    for i in range(1,len(a)):
        b[i] = 1.0/(1.0 + np.var(a)/a[i-1]**2)
    return b

def MISE(a):
    """
    Compute the Mean Integrated Square Error (MISE) defined in
    Weinberg96 Eq.7.

    Input:
    ------
    a: Coefficients.

    Output:
    -------
    D: The sum of all the coeffice
    """
    D = np.zeros(len(a))
    b = smooth_b(a)
    # This for is over all the coefficients
    for i in range(len(a)):
        D[i] = np.sum(b[:i+1]**2.0 * np.var(a) + (b[:i+1]-1.0)**2.0*a[:i+1]**2.0)
    return D


# ---------------------------  PCA -------------------------------
def S_matrix(a):
    a = a.flatten()
    S = np.zeros((len(a), len(a)))
    for i in range(len(a)):
        for j in range(len(a)):
            S[i][j] = a[i]*a[j]
    return S

def lambda_prime(lambdas):
    F = np.zeros(len(lambdas))
    for i in range(len(F)):
        F[i] = np.sum(lambdas.real[:i])/np.sum(lambdas.real)
    return F


def a_prime(T, a):
    a_new = np.zeros(len(a))
    for i in range(len(a_new)):
        a_new[i] = np.sum(np.dot(T.real[i],a))
    return a_new

