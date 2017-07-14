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


def load_plasma(ncells, domain_min, domain_max, injected_plasma_density, injected_plasma_ppc, plasma_min, plasma_max):
    '''

    Sets initial conditions for the plasma acceleration setup.

    '''

    comm.Barrier()

    dx = (domain_max - domain_min) / ncells

    nplasma_cells = ((plasma_max - plasma_min)/dx + 0.5).astype('l') + 1
    iplasma_min = ((plasma_min - domain_min)/dx).astype('l')*dx + domain_min

    #  Species 0 - the plasma
    plasma_weight = injected_plasma_density * dx[0] * dx[1] * dx[2] / injected_plasma_ppc

    # Fill the entire domain with particles. Only a subset of these points
    # will be selected for each species
    Z, Y, X, P = np.mgrid[0:nplasma_cells[2], 0:nplasma_cells[1], 0:nplasma_cells[0], 0:injected_plasma_ppc]
    Z = Z.flatten()
    Y = Y.flatten()
    X = X.flatten()
    P = P.flatten()

    particle_shift = (0.5 + P) / injected_plasma_ppc

    xp = (X + particle_shift)*dx[0] + iplasma_min[0]
    yp = (Y + particle_shift)*dx[1] + iplasma_min[1]
    zp = (Z + particle_shift)*dx[2] + iplasma_min[2]

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
    plasma_wp = np.full(Np, plasma_weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    electrons.add_particles(x=plasma_xp[lo:hi], y=plasma_yp[lo:hi], z=plasma_zp[lo:hi],
                            ux=plasma_uxp[lo:hi], uy=plasma_uyp[lo:hi], uz=plasma_uzp[lo:hi], 
                            w=plasma_wp[lo:hi], unique_particles=True)

    comm.Barrier()

def inject_new_plasma(ncells, domain_min, domain_max, num_shift, direction, injected_plasma_density, injected_plasma_ppc, plasma_min, plasma_max):
    '''

    This injects fresh plasma into the domain after the moving window has been upated.

    '''

    comm.Barrier()

    dx = (domain_max - domain_min) / ncells

    plasma_min_inject = plasma_min.copy()
    plasma_max_inject = plasma_max.copy()
    if (num_shift > 0):
        plasma_min_inject[direction] = domain_max[direction]
        plasma_max_inject[direction] = domain_max[direction] + num_shift*dx[direction]
    else:
        plasma_min_inject[direction] = domain_min[direction] - num_shift*dx[direction]
        plasma_max_inject[direction] = domain_min[direction]

    load_plasma(ncells, domain_min, domain_max, injected_plasma_density, injected_plasma_ppc, plasma_min_inject, plasma_max_inject)

ncells = np.array([64, 64, 480])/2
domain_min = np.array([-30.e-6, -30.e-6, -56.e-6])
domain_max = np.array([ 30.e-6,  30.e-6,  12.e-6])
dx = (domain_max - domain_min) / ncells

moving_window_velocity = np.array([0., 0., PICMI.clight])

injected_plasma_density = 1.e23
injected_plasma_ppc = 4
plasma_min = np.array([-20.e-6, -20.e-6,  0.0e-6])
plasma_max = np.array([ 20.e-6,  20.e-6,  domain_max[2]])

grid = PICMI.Grid(nx=ncells[0], ny=ncells[1], nz=ncells[2],
                  xmin=domain_min[0], xmax=domain_max[0], ymin=domain_min[1], ymax=domain_max[1], zmin=domain_min[2], zmax=domain_max[2],
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='open', bczmax='open',
                  moving_window_velocity = moving_window_velocity,
                  max_grid_size=32, max_level=0, coord_sys=0)

solver = PICMI.EM_solver(current_deposition = 3,
                         charge_deposition = 0,
                         field_gathering = 0,
                         particle_pusher = 0)

laser = PICMI.Gaussian_laser(antenna_z0 = 9.e-6,  # This point is on the laser plane
                             antenna_zvec = 1.,  # The plane normal direction
                             pol_angle = PICMI.pi/2.,  # The main polarization vector
                             E0 = 16.e12,  # Maximum amplitude of the laser field (in V/m)
                             waist = 5.e-6,  # The waist of the laser (in meters)
                             duration = 15.e-15,  # The duration of the laser (in seconds)
                             t_peak = 30.e-15,  # The time at which the laser reaches its peak (in seconds)
                             focal_position = 109.e-6,  # Focal position (m)
                             wavelength = 0.8e-6,  # The wavelength of the laser (in meters)
                             em_solver = solver)

electrons = PICMI.Species(type=PICMI.Electron, name='electrons')

# Maximum number of time steps
max_step = 1000

sim = PICMI.Simulation(plot_int = 200,
                       verbose = 1,
                       cfl = 1.0)

load_plasma(ncells, domain_min, domain_max, injected_plasma_density, injected_plasma_ppc, plasma_min, plasma_max)

direction = np.argmax(abs(moving_window_velocity))
old_x = grid.getmins()[direction]
new_x = old_x
for i in range(1, max_step + 1):

    # check whether the moving window has updated
    num_shift = int( 0.5 + (new_x - old_x) / dx[direction])
    if (num_shift != 0):
        inject_new_plasma(ncells, domain_min, domain_max, num_shift, direction, injected_plasma_density, injected_plasma_ppc, plasma_min, plasma_max)
        domain_min[direction] += num_shift*dx[direction]
        domain_max[direction] += num_shift*dx[direction]

    sim.step(1)

    old_x = new_x
    new_x = grid.getmins()[direction]

sim.finalize()
