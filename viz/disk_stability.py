import numpy as np
import yt
import sys
from pygadgetreader import readsnap
import octopus
from yt.units import kpc, Msun, parsec
"""
input: snap_name
return: X-Y scatter plot with color map corresponding to the particles
mass.

"""

# Reading data
snap_name = sys.argv[1]
out_name = sys.argv[2]

disk_pos = readsnap('../data/MW_models/' + snap_name, 'pos', 'disk')
disk_vel = readsnap('../data/MW_models/' + snap_name, 'vel', 'disk')
disk_mass = readsnap('../data/MW_models/' + snap_name, 'mass', 'disk')
disk_pot = readsnap('../data/MW_models/' + snap_name, 'pot', 'disk')

# computing COM with Octopus.
pos_cm, vel_cm = octopus.CM(disk_pos, disk_vel)
print('COM:' , pos_cm, vel_cm)

# Truncating disk particles
r = ((disk_pos[:,0] - disk_cm[0])**2 + (halo_disk[:,1]-pos_cm[1])**2 + (disk_pos[:,2]-pos_cm[2])**2)**0.5

index_cut = np.where(r<25)[0]


# transforming to disk COM positions
disk_pos = disk_pos[index_cut]
disk_vel = disk_vel[index_cut]
disk_pot = disk_pot[index_cut]
disk_mass = disk_mass[index_cut]


ppx = disk_pos[:,0] - pos_cm[0]
ppy = disk_pos[:,1] - pos_cm[1]
ppz = disk_pos[:,2] - pos_cm[2]

data = {'particle_position_x': ppx,\
        'particle_position_y': ppy,\
        'particle_position_z': ppz,\
        'particle_velocity_x': disk_vel[:,0] - vel_cm[0],\
        'particle_velocity_y': disk_vel[:,1] - vel_cm[1],\
        'particle_velocity_z': disk_vel[:,2] - vel_cm[2],\
        'particle_mass': disk_mass,\
        'particle_pot' : disk_pot}

# Defining box dimensions.
bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])

#ds = yt.load_particles(data, length_unit=kpc, mass_unit=1e10*Msun, n_ref=32, bbox=bbox)
ds = yt.load_particles(data, length_unit=kpc, mass_unit=1e10*Msun, bbox=bbox)
#ds.index
ad = ds.all_data()

print(ad)

#print(ds.derived_field_list)


cic_density = ad["deposit", "all_cic"]
nn_density = ad["deposit", "all_density"]
nn_deposited_mass = ad["deposit", "all_mass"]
particle_count_per_cell = ad["deposit", "all_count"]

#p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass')

p = yt.ParticleProjectionPlot(ds, 2, ['particle_mass'], depth=(5, 'kpc'))
#slc = yt.SlicePlot(ds, 2, ('deposit', 'all_cic'))
#slc.set_axes_unit('kpc')
#slc.set_unit('particle_mass', 'Msun')
#lc.set_xlim(-50, 50)
#slc.set_ylim(-50, 50)
#slc.zoom(10)
#slc.save(out_name)


p.set_width((30, 'kpc'))
p.set_axes_unit('kpc')
p.save(out_name)
