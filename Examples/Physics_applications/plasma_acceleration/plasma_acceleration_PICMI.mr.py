import numpy as np
from pywarpx import picmi
import pywarpx
#from warp import picmi

nx = 64
ny = 64
nz = 64

xmin = -200.e-6
xmax = +200.e-6
ymin = -200.e-6
ymax = +200.e-6
zmin = -200.e-6
zmax = +200.e-6

moving_window_velocity = [0., 0., picmi.c]

number_per_cell_each_dim = [4, 4, 4]

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                             moving_window_velocity = moving_window_velocity,
                             #refined_regions = [[1, [-25e-6, -25e-6, -200.e-6], [25e-6, 25e-6, 200.e-6]]],  # as argument
                             warpx_max_grid_size=128, warpx_blocking_factor=16)

# --- As a seperate function call (instead of refined_regions argument)
grid.add_refined_region(level = 1,
                        lo = [-25e-6, -25e-6, -200.e-6],
                        hi = [25e-6, 25e-6, 200.e-6])

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1,
                                     warpx_do_pml = 1,
                                     warpx_pml_ncell = 10)

beam_distribution = picmi.UniformDistribution(density = 1.e23,
                                              lower_bound = [-20.e-6, -20.e-6, -150.e-6],
                                              upper_bound = [+20.e-6, +20.e-6, -100.e-6],
                                              directed_velocity = [0., 0., 1.e9])

plasma_distribution = picmi.UniformDistribution(density = 1.e22,
                                                lower_bound = [-200.e-6, -200.e-6, 0.],
                                                upper_bound = [+200.e-6, +200.e-6, None],
                                                fill_in = True)

beam = picmi.Species(particle_type='electron', name='beam', initial_distribution=beam_distribution)
plasma = picmi.Species(particle_type='electron', name='plasma', initial_distribution=plasma_distribution)

sim = picmi.Simulation(solver = solver,
                       max_steps = 1000,
                       verbose = 1,
                       warpx_plot_int = 2,
                       warpx_current_deposition_algo = 3,
                       warpx_charge_deposition_algo = 0,
                       warpx_field_gathering_algo = 0,
                       warpx_particle_pusher_algo = 0)

sim.add_species(beam, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))
sim.add_species(plasma, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name = 'inputs_from_PICMI.mr')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

