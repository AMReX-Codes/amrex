import numpy as np
from pywarpx import picmi

constants = picmi.constants

nx = 64
ny = 64

xmin = -20.e-6
ymin = -20.e-6
xmax = +20.e-6
ymax = +20.e-6

uniform_plasma = picmi.UniformDistribution(density = 1.e25,
                                           upper_bound = [0., None, None],
                                           directed_velocity = [0.1*constants.c, 0., 0.])

electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=uniform_plasma)

grid = picmi.Cartesian2DGrid(number_of_cells = [nx, ny],
                             lower_bound = [xmin, ymin],
                             upper_bound = [xmax, ymax],
                             lower_boundary_conditions = ['periodic', 'periodic'],
                             upper_boundary_conditions = ['periodic', 'periodic'],
                             moving_window_velocity = [0., 0., 0.],
                             warpx_max_grid_size=32)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

sim = picmi.Simulation(solver = solver,
                       max_steps = 40,
                       verbose = 1,
                       warpx_plot_int = 1,
                       warpx_current_deposition_algo = 'direct')

sim.add_species(electrons, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2,2], grid=grid))

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name='inputs2d_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

