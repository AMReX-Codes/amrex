import numpy as np
from pywarpx import PICMI
#from warp import PICMI
if PICMI.codename == 'WarpX':
    from pywarpx import PGroups

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


def set_initial_conditions(ncells, domain_min, domain_max):
    '''

    Sets initial conditions for the plasma acceleration setup.


    '''

    comm.Barrier()

    num_ppc = 4
    c = 299792458.0
    dx = (domain_max - domain_min) / ncells

    #  Species 0 - the beam
    beam_min = np.array([ -20e-6,  -20e-6, -150e-6])
    beam_max = np.array([  20e-6,   20e-6, -100e-6])
    beam_density = 1.e23
    beam_gamma = 1.e9
    beam_uz = np.sqrt(beam_gamma**2 - 1)*c
    beam_weight = beam_density * dx[0] * dx[1] * dx[2] / num_ppc

    #  Species 1 - the plasma
    plasma_min = np.array([-200e-6, -200e-6,  0.0e-6])
    plasma_max = np.array([ 200e-6,  200e-6,  200e-6])
    plasma_density = 1.e22
    plasma_weight = plasma_density * dx[0] * dx[1] * dx[2] / num_ppc

    # Fill the entire domain with particles. Only a subset of these points
    # will be selected for each species
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
    beam_wp = np.full(Np, beam_weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    beam.add_particles(x=beam_xp[lo:hi], y=beam_yp[lo:hi], z=beam_zp[lo:hi],
                       ux=beam_uxp[lo:hi], uy=beam_uyp[lo:hi], uz=beam_uzp[lo:hi], 
                       w=beam_wp[lo:hi], unique_particles=1)

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

    plasma.add_particles(x=plasma_xp[lo:hi], y=plasma_yp[lo:hi], z=plasma_zp[lo:hi],
                         ux=plasma_uxp[lo:hi], uy=plasma_uyp[lo:hi], uz=plasma_uzp[lo:hi], 
                         w=plasma_wp[lo:hi], unique_particles=1)

    comm.Barrier()


def inject_plasma(num_shift, direction):
    '''

    This injects fresh plasma into the domain after the moving window has been upated.

    '''

    comm.Barrier()

    num_ppc = 4
    density = 1.e22
    weight = density * dx[0]*dx[1]*dx[2] / num_ppc

    shift_box = np.array(ncells)
    shift_box[direction] = abs(num_shift)

    Z, Y, X, P = np.mgrid[0:shift_box[2], 0:shift_box[1], 0:shift_box[0], 0:num_ppc]
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()
    P = P.flatten()
    particle_shift = (0.5 + P) / num_ppc

    pos = [[],[],[]]
    
    if (num_shift > 0):
        pos[0] = (X + particle_shift)*dx[0] + domain_min[0]
        pos[1] = (Y + particle_shift)*dx[1] + domain_min[1]
        pos[2] = (Z + particle_shift)*dx[2] + domain_min[2]
        pos[direction] -= domain_min[direction]
        pos[direction] += domain_max[direction]
    
    if (num_shift < 0):
        pos[0] = (X + particle_shift)*dx[0] + domain_min[0]
        pos[1] = (Y + particle_shift)*dx[1] + domain_min[1]
        pos[2] = (Z + particle_shift)*dx[2] + domain_min[2]
        pos[direction] += num_shift*dx[direction]

    xp = pos[0]
    yp = pos[1]
    zp = pos[2]

    uxp = np.zeros_like(xp)
    uyp = np.zeros_like(xp)
    uzp = np.zeros_like(xp)

    Np = xp.shape[0]
    wp = np.full(Np, weight)

    lo, hi = get_parallel_indices(Np, comm.rank, comm.size)

    plasma.add_particles(x=xp[lo:hi], y=yp[lo:hi], z=zp[lo:hi],
                         ux=uxp[lo:hi], uy=uyp[lo:hi], uz=uzp[lo:hi], 
                         w=wp[lo:hi], unique_particles=1)

    comm.Barrier()

ncells = np.array([64, 64, 64])
domain_min = np.array([-200e-6, -200e-6, -200e-6])
domain_max = np.array([ 200e-6,  200e-6,  200e-6])
dx = (domain_max - domain_min) / ncells

moving_window_velocity = np.array([0., 0., PICMI.clight])

grid = PICMI.Grid(nx=ncells[0], ny=ncells[1], nz=ncells[2],
                  xmin=domain_min[0], xmax=domain_max[0], ymin=domain_min[1], ymax=domain_max[1], zmin=domain_min[2], zmax=domain_max[2],
                  bcxmin='periodic', bcxmax='periodic', bcymin='periodic', bcymax='periodic', bczmin='open', bczmax='open',
                  moving_window_velocity = moving_window_velocity,
                  max_grid_size=32, max_level=0, coord_sys=0)


solver = PICMI.EM_solver(current_deposition = 3,
                         charge_deposition = 0,
                         field_gathering = 0,
                         particle_pusher = 0)

beam = PICMI.Species(type=PICMI.Electron, name='beam')
plasma = PICMI.Species(type=PICMI.Electron, name='plasma')


# Maximum number of time steps
max_step = 160

sim = PICMI.Simulation(plot_int = 2,
                       verbose = 1,
                       cfl = 1.0)

set_initial_conditions(ncells, domain_min, domain_max)

direction = np.argmax(abs(moving_window_velocity))
old_x = grid.getmins()[direction]
new_x = old_x
for i in range(1, max_step + 1):

    # check whether the moving window has updated
    num_shift = int( 0.5 + (new_x - old_x) / dx[direction])
    if (num_shift != 0):
        inject_plasma(num_shift, direction)
        domain_min[direction] += num_shift*dx[direction]
        domain_max[direction] += num_shift*dx[direction]

    sim.step(1)

    old_x = new_x
    new_x = grid.getmins()[direction]

sim.finalize()
