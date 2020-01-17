# --- Simple example of Langmuir oscillations in a uniform plasma

from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

plasma_density = 1.e25
plasma_xmin = 0.
plasma_x_velocity = 0.1*constants.c

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 40
diagnostic_interval = 10

# --- Grid
nx = 64
ny = 64
nz = 64

xmin = -20.e-6
ymin = -20.e-6
zmin = -20.e-6
xmax = +20.e-6
ymax = +20.e-6
zmax = +20.e-6

number_per_cell_each_dim = [2,2,2]

##########################
# physics components
##########################

uniform_plasma = picmi.UniformDistribution(density = 1.e25,
                                           upper_bound = [0., None, None],
                                           directed_velocity = [0.1*constants.c, 0., 0.])

electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=uniform_plasma)

##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'periodic'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'periodic'],
                             moving_window_velocity = [0., 0., 0.],
                             warpx_max_grid_size = 32)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(grid = grid,
                                    period = diagnostic_interval,
                                    data_list = ['Ex', 'Jx'])

part_diag1 = picmi.ParticleDiagnostic(period = diagnostic_interval,
                                      species = [electrons],
                                      data_list = ['weighting', 'ux', 'Ex'])

##########################
# simulation setup
##########################

sim = picmi.Simulation(solver = solver,
                       max_steps = max_steps,
                       verbose = 1,
                       warpx_current_deposition_algo = 'direct')

sim.add_species(electrons,
                layout = picmi.GriddedLayout(n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid))

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
#sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()

