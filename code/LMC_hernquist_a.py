"""
This code find the scale length 'a'
of the Hernquist profile in order to
match the a desired enclosed mass M_ideal.

Input:
     M: total mass
     M_ideal: Desired enclosed mass

Output:
     a: Hernquist scale length.
"""
import numpy as np
from soda import profiles as pfls
import sys

#M = float(sys.argv[1])
#M_ideal = float(sys.argv[2])

def a_value(M, M_ideal):
    # possible a candidates.
    a_values = np.linspace(3, 40, 2000)
    best_a = []
    # r_obs is where the rotation curve peaks.
    r_obs = 8.7 # kpc

    for i in a_values:
        Mass_enclosed = pfls.mass_hernquist(i, r_obs, M) #Enclosed mass within r_obs
        if (abs(Mass_enclosed-M_ideal) < 5E8):
            best_a = i
    return best_a
