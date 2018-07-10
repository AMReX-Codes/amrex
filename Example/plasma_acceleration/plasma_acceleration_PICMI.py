import numpy as np
from pywarpx import PICMI
#from warp import PICMI

nx = 64
ny = 64
nz = 64

xmin = -200.e-6
xmax = +200.e-6
ymin = -200.e-6
ymax = +200.e-6
zmin = -200.e-6
zmax = +200.e-6

moving_window_velocity = [0., 0., PICMI.c]

number_per_cell_each_dim = [2, 2, 1]

grid = PICMI.Grid(nx=nx, ny=ny, nz=nz, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='open', bczmax='open',
                  moving_window_velocity = moving_window_velocity,
                  max_grid_size=32, max_level=0, coord_sys=0)

solver = PICMI.EM_solver(current_deposition_algo = 3,
                         charge_deposition_algo = 0,
                         field_gathering_algo = 0,
                         particle_pusher_algo = 0)

beam = PICMI.Species(type='electron', name='beam')
plasma = PICMI.Species(type='electron', name='plasma')

beam_distribution = PICMI.Plasma(species = beam,
                                 density = 1.e23,
                                 xmin = -20.e-6, xmax = +20.e-6,
                                 ymin = -20.e-6, ymax = +20.e-6,
                                 zmin = -150.e-6, zmax = -100.e-6,
                                 vzmean = 1.e9,
                                 number_per_cell_each_dim = number_per_cell_each_dim)

plasma_distribution = PICMI.Plasma(species = plasma,
                                   density = 1.e22,
                                   xmin = -200.e-6, xmax = +200.e-6,
                                   ymin = -200.e-6, ymax = +200.e-6,
                                   zmin = 0.,
                                   number_per_cell_each_dim = number_per_cell_each_dim,
                                   fill_in = True)

sim = PICMI.Simulation(plot_int = 2,
                       verbose = 1,
                       cfl = 1.0,
                       max_step = 1000)

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_inputs(inputs_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

