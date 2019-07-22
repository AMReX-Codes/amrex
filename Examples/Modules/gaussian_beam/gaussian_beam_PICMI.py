import numpy as np
from pywarpx import picmi
#from warp import picmi

nx = 32
ny = 32
nz = 32

xmin = -2.
xmax = +2.
ymin = -2.
ymax = +2.
zmin = -2.
zmax = +2.

number_sim_particles = 32768
total_charge = 8.010883097437485e-07

beam_rms_size = 0.25
electron_beam_divergence = -0.04*picmi.c

em_order = 3

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                             warpx_max_grid_size=16)

solver = picmi.ElectromagneticSolver(grid = grid,
                                     cfl = 1.,
                                     stencil_order=[em_order,em_order,em_order])

electron_beam = picmi.GaussianBunchDistribution(n_physical_particles = total_charge/picmi.q_e,
                                                rms_bunch_size = [beam_rms_size, beam_rms_size, beam_rms_size],
                                                velocity_divergence = [electron_beam_divergence, electron_beam_divergence, electron_beam_divergence])

proton_beam = picmi.GaussianBunchDistribution(n_physical_particles = total_charge/picmi.q_e,
                                              rms_bunch_size = [beam_rms_size, beam_rms_size, beam_rms_size])

electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=electron_beam)
protons = picmi.Species(particle_type='proton', name='protons', initial_distribution=proton_beam)

sim = picmi.Simulation(solver = solver,
                       max_steps = 1000,
                       verbose = 1,
                       warpx_plot_int = 8,
                       warpx_current_deposition_algo = 0,
                       warpx_charge_deposition_algo = 0,
                       warpx_field_gathering_algo = 0,
                       warpx_particle_pusher_algo = 0)

sim.add_species(electrons, layout=picmi.PseudoRandomLayout(n_macroparticles=number_sim_particles))
sim.add_species(protons, layout=picmi.PseudoRandomLayout(n_macroparticles=number_sim_particles))

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

