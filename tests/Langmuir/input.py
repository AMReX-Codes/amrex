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
amr.n_cell =   "%d  %d  %d"%(nx, ny, nz)

# Maximum allowable size of each subdomain in the problem domain; 
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 64

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

amr.plot_int = 1   # How often to write plotfiles.  "<= 0" means no plotfiles.

# Geometry
geometry.coord_sys   = 0                  # 0: Cartesian
geometry.is_periodic = "1     1     1"      # Is periodic?  
geometry.prob_lo     = "%7.0e   %7.0e   %7.0e"%(xmin, ymin, zmin)    # physical domain
geometry.prob_hi     = "%7.0e   %7.0e   %7.0e"%(xmax, ymax, zmax)

# Verbosity
warpx.verbose = 1

# Algorithms
algo.current_deposition = 3
algo.charge_deposition = 0
algo.field_gathering = 1
algo.particle_pusher = 0

# CFL
warpx.cfl = 1.0

# --- Initialize the simulation
boxlib = BoxLib()
boxlib.init()
warpx.init()

# --- Set initial conditions for the Languir wave setup and pass to WarpX.

n_e = 1.0e25
num_particles_per_cell = 10
ux = 0.01
uy = 0.0
uz = 0.0

weight = n_e * dx * dy * dz / num_particles_per_cell
gamma = 1.0 / np.sqrt(1.0 - ux**2 - uy**2 - uz**2)
c = 299792458.0
ux *= gamma*c
uy *= gamma*c
uz *= gamma*c

# --- This creates particles only the left half of the domain
Z, Y, X, P = np.mgrid[0:64, 0:64, 0:32, 0:num_particles_per_cell]
Z = Z.flatten()
Y = Y.flatten()
X = X.flatten()
P = P.flatten()

particle_shift = (0.5 + P) / num_particles_per_cell

# -- particle positions
xp = (X + particle_shift)*dx + xmin
yp = (Y + particle_shift)*dy + ymin
zp = (Z + particle_shift)*dz + zmin

# --- velocities
uxp = ux * np.ones_like(xp)
uyp = uy * np.ones_like(xp)
uzp = uz * np.ones_like(xp)

# --- and weights
wp = np.full((xp.shape[0], 1), weight)

warpxC.addNParticles(0, xp, yp, zp, uxp, uyp, uzp, wp, 1)

warpx.evolve(max_step)
warpx.finalize()

boxlib.finalize()
