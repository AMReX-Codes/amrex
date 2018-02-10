import numpy as np
from pywarpx import PICMI

nx = 64
ny = 64
nz = 64

xmin = -20.e-6
ymin = -20.e-6
zmin = -20.e-6
xmax = +20.e-6
ymax = +20.e-6
zmax = +20.e-6

dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny
dz = (zmax - zmin)/nz

grid = PICMI.Grid(nx=nx, ny=ny, nz=nz,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='periodic', bczmax='periodic',
                  max_grid_size=32, max_level=0, coord_sys=0)

solver = PICMI.EM_solver(current_deposition_algo = 3,
                         charge_deposition_algo = 0,
                         field_gathering_algo = 0,
                         particle_pusher_algo = 0)

electrons = PICMI.Species(type=PICMI.Electron, name='electrons', lvariableweights=True)

plasma = PICMI.Plasma(species=electrons, density=1.e25, xmax=0., number_per_cell_each_dim=[2,2,2])

sim = PICMI.Simulation(plot_int = 1,
                       verbose = 1,
                       cfl = 1.0,
                       max_step = 40)

sim.initialize(inputs_name='inputs_from_PICMI')

