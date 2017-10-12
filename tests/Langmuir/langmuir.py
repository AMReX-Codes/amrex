import numpy as np
from pywarpx import *

from mpi4py import MPI
comm = MPI.COMM_WORLD

def get_parallel_indices(Np, rank, size):
    '''

    This decomposes the arrays of particle attributes into subarrays
    for parallel processing.

    '''
    Navg = Np / size
    Nleft = Np - Navg * size
    if (rank < Nleft):
        lo = rank*(Navg + 1)
        hi = lo + Navg + 1
    else:
        lo = rank * Navg + Nleft
        hi = lo + Navg
    return lo, hi


def set_initial_conditions():
    '''

    Sets initial conditions for the Langmuir wave setup.


    '''

    comm.Barrier()

    n_e = 1.0e25
    num_particles_per_cell = 10
    ux = 0.01
    uy = 0.0
    uz = 0.0

    weight = n_e * dx * dy * dz / num_particles_per_cell
    c = 299792458.0
    ux *= c
    uy *= c
    uz *= c

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
    Np = xp.shape[0]
    wp = np.full((Np, 1), weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    add_particles(0, xp[lo:hi], yp[lo:hi], zp[lo:hi],
                  uxp[lo:hi], uyp[lo:hi], uzp[lo:hi], 
                  wp[lo:hi], 1)

    comm.Barrier()

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
amr.max_grid_size = 32

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

amr.plot_int = 10   # How often to write plotfiles.  "<= 0" means no plotfiles.

# Geometry
geometry.coord_sys   = 0                  # 0: Cartesian
geometry.is_periodic = "1     1     1"      # Is periodic?  
geometry.prob_lo     = "%7.0e   %7.0e   %7.0e"%(xmin, ymin, zmin)    # physical domain
geometry.prob_hi     = "%7.0e   %7.0e   %7.0e"%(xmax, ymax, zmax)

warpx.serialize_ics = 1

# Verbosity
warpx.verbose = 1

particles.nspecies = 1
particles.species_names = "electrons"

# Algorithms
algo.current_deposition = 3
algo.charge_deposition = 0
algo.field_gathering = 1
algo.particle_pusher = 0

# CFL
warpx.cfl = 1.0

# --- Initialize the simulation
amrex = AMReX()
amrex.init()
warpx.init()

set_initial_conditions()

warpx.evolve(max_step)
warpx.finalize()

amrex.finalize()
