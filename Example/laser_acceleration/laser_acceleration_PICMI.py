import numpy as np
from pywarpx import picmi
#from warp import picmi

nx = 64
ny = 64
nz = 480

xmin = -30.e-6
xmax = +30.e-6
ymin = -30.e-6
ymax = +30.e-6
zmin = -56.e-6
zmax = +12.e-6

moving_window_velocity = [0., 0., picmi.c]

injected_plasma_density = 1.e23
number_per_cell_each_dim = [2, 2, 1]
plasma_min = [-20.e-6, -20.e-6,  0.0e-6]
plasma_max = [ 20.e-6,  20.e-6,  zmax]

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                             moving_window_velocity = moving_window_velocity,
                             max_grid_size=32, max_level=0, coord_sys=0)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

t_peak = 30.e-15  # The time at which the laser reaches its peak at the antenna (in seconds)
focal_distance = 100.e-6  # Focal distance from the antenna (in meters)
antenna_z0 = 9.e-6  # This point is on the laser plane
laser = picmi.GaussianLaser(wavelength = 0.8e-6,  # The wavelength of the laser (in meters)
                            waist = 5.e-6,  # The waist of the laser (in meters)
                            duration = 15.e-15,  # The duration of the laser (in seconds)
                            polarization_angle = np.pi/2.,  # The main polarization vector
                            focal_position = [0., 0., focal_distance + antenna_z0],  # Focal position (m)
                            E0 = 16.e12,  # Maximum amplitude of the laser field (in V/m)
                            centroid_position = [0., 0., antenna_z0 - picmi.c*t_peak], # Position of the laser centroid in Z at time 0
                            propagation_direction = [0,0,1])

laser_antenna = picmi.LaserAntenna(position = [0., 0., antenna_z0],  # This point is on the laser plane
                                   normal_vector = [0., 0., 1.])  # The plane normal direction

uniform_plasma = picmi.UniformDistribution(density = injected_plasma_density,
                                           lower_bound = plasma_min,
                                           upper_bound = plasma_max,
                                           fill_in = True)

electrons = picmi.Species(particle_type = 'electron',
                          name = 'electrons',
                          initial_distribution = uniform_plasma)

sim = picmi.Simulation(solver = solver,
                       plot_int = 100,
                       verbose = 1,
                       cfl = 1.0,
                       max_steps = 1000,
                       current_deposition_algo = 3,
                       charge_deposition_algo = 0,
                       field_gathering_algo = 0,
                       particle_pusher_algo = 0)

sim.add_species(electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))

sim.add_laser(laser, injection_method=laser_antenna)

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

