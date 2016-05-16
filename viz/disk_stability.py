import numpy as np
import yt
import sys
from pygadgetreader import *

snap_name = sys.argv[1]

bbox = [[-89000, 89000],
        [-89000, 89000],
        [-89000, 89000]]


ds = yt.load('../data/MW_models/pm/UHR/' + snap_name, bounding_box=bbox)
ds.index
ad =  ds.all_data()

p = yt.ParticlePlot(ds, ('PartType2','particle_position_x'),
('PartType2','particle_position_y'), ('PartType2','particle_mass'))
p.set_unit('particle_mass', 'Msun')
p.set_xlim(-20, 20)
p.set_ylim(-20, 20)
p.save('MW2_disk.png')
