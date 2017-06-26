import numpy as np
from pywarpx import PICMI
#from warp import PICMI

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
    wp = np.full(Np, weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    electrons.add_particles(x=xp[lo:hi], y=yp[lo:hi], z=zp[lo:hi],
                            ux=uxp[lo:hi], uy=uyp[lo:hi], uz=uzp[lo:hi], 
                            w=wp[lo:hi], unique_particles=1)

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

grid = PICMI.Grid(nx=nx, ny=ny, nz=nz,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='periodic', bczmax='periodic',
                  max_grid_size=32, max_level=0, coord_sys=0)

solver = PICMI.EM_solver(current_deposition = 3,
                         charge_deposition = 0,
                         field_gathering = 1,
                         particle_pusher = 0)

electrons = PICMI.Species(type=PICMI.Electron, name='electrons', lvariableweights=True)

sim = PICMI.Simulation(plot_int = 10,
                       verbose = 1,
                       cfl = 1.0)

set_initial_conditions()

# Maximum number of time steps
max_step = 40

sim.step(max_step)

sim.finalize()
