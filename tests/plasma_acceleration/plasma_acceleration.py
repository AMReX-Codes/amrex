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

    Sets initial conditions for the plasma acceleration setup.


    '''

    comm.Barrier()

    ncells = np.array([nx, ny, nz])
    num_ppc = 4

    domain_min = np.array([xmin, ymin, zmin])
    domain_max = np.array([xmax, ymax, zmax])

    beam_min = np.array([ -20e-6,  -20e-6, -150e-6])
    beam_max = np.array([  20e-6,   20e-6, -100e-6])
    
    plasma_min = np.array([-200e-6, -200e-6,  0.0e-6])
    plasma_max = np.array([ 200e-6,  200e-6,  200e-6])

    dx = (domain_max - domain_min) / ncells

    c = 299792458.0

    #  Species 0 - the beam
    beam_density = 1.e23
    beam_gamma = 1.e9
    beam_uz = np.sqrt(beam_gamma**2 -1)*c
    beam_weight = beam_density * dx[0] * dx[1] * dx[2] / num_ppc

    #  Species 1 - the plasma
    plasma_density = 1.e22
    plasma_weight = plasma_density * dx[0] * dx[1] * dx[2] / num_ppc

    Z, Y, X, P = np.mgrid[0:ncells[0], 0:ncells[1], 0:ncells[2], 0:num_ppc]
    Z = Z.flatten()
    Y = Y.flatten()
    X = X.flatten()
    P = P.flatten()

    particle_shift = (0.5 + P) / num_ppc

    xp = (X + particle_shift)*dx[0] + domain_min[0]
    yp = (Y + particle_shift)*dx[1] + domain_min[1]
    zp = (Z + particle_shift)*dx[2] + domain_min[2]

    # Do the beam species first
    beam_locs = np.logical_and(xp >= beam_min[0], xp < beam_max[0])
    beam_locs = np.logical_and(beam_locs, yp >= beam_min[1])
    beam_locs = np.logical_and(beam_locs, zp >= beam_min[2])
    beam_locs = np.logical_and(beam_locs, yp  <  beam_max[1])
    beam_locs = np.logical_and(beam_locs, zp  <  beam_max[2])

    beam_xp = xp[beam_locs]
    beam_yp = yp[beam_locs]
    beam_zp = zp[beam_locs]

    beam_uxp = np.zeros_like(beam_xp)
    beam_uyp = np.zeros_like(beam_xp)
    beam_uzp = beam_uz * np.ones_like(beam_xp)

    # --- and weights
    Np = beam_xp.shape[0]
    beam_wp = np.full((Np, 1), beam_weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    warpxC.addNParticles(0, beam_xp[lo:hi], beam_yp[lo:hi], beam_zp[lo:hi],
                         beam_uxp[lo:hi], beam_uyp[lo:hi], beam_uzp[lo:hi], 
                         beam_wp[lo:hi], 1)

    # now do the plasma species
    plasma_locs = np.logical_and(xp >= plasma_min[0], xp < plasma_max[0])
    plasma_locs = np.logical_and(plasma_locs, yp >= plasma_min[1])
    plasma_locs = np.logical_and(plasma_locs, zp >= plasma_min[2])
    plasma_locs = np.logical_and(plasma_locs, yp  <  plasma_max[1])
    plasma_locs = np.logical_and(plasma_locs, zp  <  plasma_max[2])

    plasma_xp = xp[plasma_locs]
    plasma_yp = yp[plasma_locs]
    plasma_zp = zp[plasma_locs]

    plasma_uxp = np.zeros_like(plasma_xp)
    plasma_uyp = np.zeros_like(plasma_xp)
    plasma_uzp = np.zeros_like(plasma_xp)

    # --- and weights
    Np = plasma_xp.shape[0]
    plasma_wp = np.full((Np, 1), plasma_weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    warpxC.addNParticles(1, plasma_xp[lo:hi], plasma_yp[lo:hi], plasma_zp[lo:hi],
                         plasma_uxp[lo:hi], plasma_uyp[lo:hi], plasma_uzp[lo:hi], 
                         plasma_wp[lo:hi], 1)

    comm.Barrier()

nx = 64
ny = 64
nz = 64

xmin = -200.e-6
ymin = -200.e-6
zmin = -200.e-6
xmax = +200.e-6
ymax = +200.e-6
zmax = +200.e-6

dx = (xmax - xmin)/nx
dy = (ymax - ymin)/ny
dz = (zmax - zmin)/nz

# Maximum number of time steps
max_step = 60

# number of grid points
amr.n_cell =   "%d  %d  %d" % (nx, ny, nz)

# Maximum allowable size of each subdomain in the problem domain; 
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 32

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 0

amr.plot_int = 2   # How often to write plotfiles.  "<= 0" means no plotfiles.

# Geometry
geometry.coord_sys   = 0                  # 0: Cartesian
geometry.is_periodic = "1     1     0"      # Is periodic?  
geometry.prob_lo     = "%7.0e   %7.0e   %7.0e" % (xmin, ymin, zmin)    # physical domain
geometry.prob_hi     = "%7.0e   %7.0e   %7.0e" % (xmax, ymax, zmax)

# Algorithms
algo.current_deposition = 3
algo.charge_deposition = 0
algo.field_gathering = 0
algo.particle_pusher = 0

# Particles
particles.nspecies = 2

warpx.verbose = 1
warpx.cfl = 1.0
warpx.do_moving_window = 1
warpx.moving_window_dir = 2
warpx.moving_window_v = 1.0  # in units of the speed of light

warpx.do_plasma_injection = 1
warpx.num_injected_species = 1
warpx.injected_plasma_species = 1
warpx.injected_plasma_density = 1e22
warpx.injected_plasma_ppc = 4

# --- Initialize the simulation
boxlib = BoxLib()
boxlib.init()
warpx.init()

set_initial_conditions()

for i in range(max_step + 1):
    warpx.evolve(i)

warpx.finalize()

boxlib.finalize()
