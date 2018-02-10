import numpy as np
from pywarpx import PICMI
#from warp import PICMI

nx = 64
ny = 64
nz = 480

xmin = -30.e-6
xmax = +30.e-6
ymin = -30.e-6
ymax = +30.e-6
zmin = -56.e-6
zmax = +12.e-6

moving_window_velocity = [0., 0., PICMI.clight]

injected_plasma_density = 1.e23
number_per_cell_each_dim = [2, 2, 2]
plasma_min = [-20.e-6, -20.e-6,  0.0e-6]
plasma_max = [ 20.e-6,  20.e-6,  zmax]

grid = PICMI.Grid(nx=nx, ny=ny, nz=nz, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='open', bczmax='open',
                  moving_window_velocity = moving_window_velocity,
                  max_grid_size=32, max_level=0, coord_sys=0)

solver = PICMI.EM_solver(current_deposition_algo = 3,
                         charge_deposition_algo = 0,
                         field_gathering_algo = 0,
                         particle_pusher_algo = 0)

t_peak = 30.e-15  # The time at which the laser reaches its peak at the antenna (in seconds)
focal_distance = 100.e-6  # Focal distance from the antenna (in meters)
antenna_z0 = 9.e-6  # This point is on the laser plane
laser = PICMI.Gaussian_laser(wavelength = 0.8e-6,  # The wavelength of the laser (in meters)
                             waist = 5.e-6,  # The waist of the laser (in meters)
                             duration = 15.e-15,  # The duration of the laser (in seconds)
                             pol_angle = PICMI.pi/2.,  # The main polarization vector
                             focal_position = focal_distance + antenna_z0,  # Focal position (m)
                             E0 = 16.e12,  # Maximum amplitude of the laser field (in V/m)
                             z0 = antenna_z0 - PICMI.clight*t_peak, # Position of the laser centroid in Z at time 0
                             )

laser = PICMI.Laser_antenna(laser = laser,
                            antenna_z0 = antenna_z0,  # This point is on the laser plane
                            antenna_zvec = 1.,  # The plane normal direction
                            )

electrons = PICMI.Species(type=PICMI.Electron, name='electrons')

plasma = PICMI.Plasma(species = electrons,
                      density = injected_plasma_density,
                      xmin = plasma_min[0],
                      xmax = plasma_max[0],
                      ymin = plasma_min[1],
                      ymax = plasma_max[1],
                      zmin = plasma_min[2],
                      zmax = plasma_max[2],
                      number_per_cell_each_dim = number_per_cell_each_dim)

sim = PICMI.Simulation(plot_int = 100,
                       verbose = 1,
                       cfl = 1.0,
                       max_step = 1000)

sim.write_inputs(inputs_name = 'inputs_from_PICMI')
