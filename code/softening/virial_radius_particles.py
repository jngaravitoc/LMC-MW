import numpy as np
from pygadgetreader import *
import sys

snap_name = sys.argv[1]
rvir = np.float(sys.argv[2])

ic_snap_pos = readsnap(snap_name, 'pos', 'dm')

r = np.sqrt(ic_snap_pos[:,0]**2 + ic_snap_pos[:,1]**2.0 +
ic_snap_pos[:,2]**2.0)

index = np.where(r<rvir)[0]

print(len(ic_snap_pos[index,0]))

