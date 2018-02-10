import numpy as np
from pywarpx import *

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

# Maximum number of time steps
max_step = 40

# number of grid points
amr.n_cell =   [nx, ny, nz]

# Maximum allowable size of each subdomain in the problem domain; 
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 32

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

amr.plot_int = 1   # How often to write plotfiles.  "<= 0" means no plotfiles.

# Geometry
geometry.coord_sys   = 0                  # 0: Cartesian
geometry.is_periodic = [1, 1, 1]          # Is periodic?  
geometry.prob_lo     = [xmin, ymin, zmin]    # physical domain
geometry.prob_hi     = [xmax, ymax, zmax]

# Verbosity
warpx.verbose = 1
warpx.do_moving_window = 0
warpx.moving_window_dir = 'x'
warpx.moving_window_v = 0.0 # in units of the speed of light

# Algorithms
algo.current_deposition = 3
algo.charge_deposition = 0
algo.field_gathering = 0
algo.particle_pusher = 0

# CFL
warpx.cfl = 1.0

particles.nspecies = 1
particles.species_names = 'electrons'

electrons.charge = '-q_e'
electrons.mass = 'm_e'
electrons.injection_style = "NUniformPerCell"
electrons.num_particles_per_cell_each_dim = [2, 2, 2]

electrons.xmin = xmin
electrons.xmax = 0.e-6
electrons.ymin = ymin
electrons.ymax = ymax
electrons.zmin = zmin
electrons.zmax = zmax

electrons.profile = 'constant'
electrons.density = 1.e25  # number of electrons per m^3

electrons.momentum_distribution_type = "constant"
electrons.ux  = 0.01   # ux = gamma*beta_x

# --- Initialize the simulation
warpx.write_inputs('inputs_from_python', max_step=max_step)

