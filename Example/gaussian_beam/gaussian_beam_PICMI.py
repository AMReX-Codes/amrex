import numpy as np
from pywarpx import PICMI
#from warp import PICMI

nx = 32
ny = 32
nz = 32

xmin = -2.
xmax = +2.
ymin = -2.
ymax = +2.
zmin = -2.
zmax = +2.

number_sim_particles = 32678
total_charge = 8.010883097437485e-07

beam_rms_size = 0.25
electron_beam_divergence = -0.04

em_order = 3

grid = PICMI.Grid(nx=nx, ny=ny, nz=nz, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='open', bczmax='open',
                  max_grid_size=16, coord_sys=0)

solver = PICMI.EM_solver(current_deposition_algo = 0,
                         charge_deposition_algo = 0,
                         field_gathering_algo = 0,
                         particle_pusher_algo = 0,
                         norderx = em_order, nordery = em_order, norderz = em_order)

electrons = PICMI.Species(type='electron', name='electrons')
protons = PICMI.Species(type='electron', name='protons')

electron_beam = PICMI.GaussianBeam(electrons,
                                   number_sim_particles = number_sim_particles,
                                   number_real_particles = total_charge/PICMI.q_e,
                                   Xrms = beam_rms_size, Yrms = beam_rms_size, Zrms = beam_rms_size,
                                   UXdiv = electron_beam_divergence, UYdiv = electron_beam_divergence, UZdiv = electron_beam_divergence)

proton_beam = PICMI.GaussianBeam(protons,
                                 number_sim_particles = number_sim_particles,
                                 number_real_particles = total_charge/PICMI.q_e,
                                 Xrms = beam_rms_size, Yrms = beam_rms_size, Zrms = beam_rms_size)

sim = PICMI.Simulation(plot_int = 8,
                       verbose = 1,
                       cfl = 1.0,
                       max_step = 1000)

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_inputs(inputs_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
#sim.step()

