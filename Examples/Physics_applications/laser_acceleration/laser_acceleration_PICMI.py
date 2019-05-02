import numpy as np
from pywarpx import picmi
#from warp import picmi

##########################
# physics parameters
##########################

# --- laser

laser_a0              = 4.        # Normalized potential vector
laser_wavelength      = 8e-07     # Wavelength of the laser (in meters)
laser_waist           = 5e-06     # Waist of the laser (in meters)
laser_duration        = 15e-15    # Duration of the laser (in seconds)
laser_polarization    = np.pi/2.  # Polarization angle (in rad)
laser_injection_loc   = 9.e-6     # Position of injection (in meters, along z)
laser_focal_distance  = 100.e-6   # Focal distance from the injection (in meters)
laser_t_peak          = 30.e-15   # The time at which the laser reaches its peak 
                                  #   at the antenna injection location (in seconds)
# --- plasma

plasma_density = 1.e24 
plasma_min     = [-20.e-6, -20.e-6,  0.0e-6]
plasma_max     = [ 20.e-6,  20.e-6,  1.e-3]


##########################
# numerics parameters
##########################

# --- Nb time steps

max_steps = 1000

# --- grid

nx = 64
ny = 64
nz = 480

xmin = 1.5*plasma_min[0]
xmax = 1.5*plasma_max[0]
ymin = 1.5*plasma_min[1]
ymax = 1.5*plasma_max[1]
zmin = -56.e-6
zmax = 12.e-6

moving_window_velocity = [0., 0., picmi.c]

number_per_cell_each_dim = [2, 2, 1]

##########################
# physics components
##########################

# --- laser

laser = picmi.GaussianLaser(wavelength            = laser_wavelength,  
                            waist                 = laser_waist,                    
                            duration              = laser_duration,              
                            focal_position        = [0., 0., laser_focal_distance + laser_injection_loc],  
                            centroid_position     = [0., 0., laser_injection_loc - picmi.c*laser_t_peak], 
                            polarization_angle    = laser_polarization,  
                            propagation_direction = [0,0,1],
                            E0 = laser_a0*2.*np.pi*picmi.m_e*picmi.c**2/(picmi.q_e*laser_wavelength)) # Maximum amplitude of the laser field (in V/m)

laser_antenna = picmi.LaserAntenna(position = [0., 0., laser_injection_loc],  # This point is on the laser plane
                                   normal_vector = [0., 0., 1.])  # The plane normal direction

# --- plasma

uniform_plasma = picmi.UniformDistribution(density     = plasma_density,
                                           lower_bound = plasma_min,
                                           upper_bound = plasma_max,
                                           fill_in     = True)

electrons = picmi.Species(particle_type = 'electron',
                          name = 'electrons',
                          initial_distribution = uniform_plasma)


##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                             moving_window_velocity = moving_window_velocity,
                             warpx_max_grid_size=32)

solver = picmi.ElectromagneticSolver(grid=grid, method='CKC', cfl=1.)


##########################
# diagnostics
##########################

field_diag1 = picmi.FieldDiagnostic(grid = grid,
                                    period = 100,
                                    warpx_plot_raw_fields = 1,
                                    warpx_plot_raw_fields_guards = 1,
                                    warpx_plot_finepatch = 1,
                                    warpx_plot_crsepatch = 1)

part_diag1 = picmi.ParticleDiagnostic(period = 100,
                                      species = [electrons])

##########################
# simulation setup
##########################

sim = picmi.Simulation(solver = solver,
                       max_steps = max_steps,
                       verbose = 1,
                       cfl = 1.0,
                       warpx_current_deposition_algo = 3,
                       warpx_charge_deposition_algo = 0,
                       warpx_field_gathering_algo = 0,
                       warpx_particle_pusher_algo = 0)

sim.add_species(electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))

sim.add_laser(laser, injection_method=laser_antenna)

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)

##########################
# simulation run
##########################

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step(max_steps)

