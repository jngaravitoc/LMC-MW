import numpy as np 
from profiles import mass_hernquist as mh
import sys

M = float(sys.argv[1])
M_ideal = 1.3E10

def a_value(M):
    a_values = np.linspace(3, 30, 1000)
    best_a = []

    for i in a_values:
        Mass_enclosed = mh(i, 9, M)
        Mass_enclosed = Mass_enclosed.value
        if (abs(Mass_enclosed-M_ideal) < 5E8):
            best_a = i
         
    return best_a

print a_value(M)
