import matplotlib
matplotlib.use('Agg')
from mpi4py import MPI
import glob, read_raw_data
import numpy as np
import scipy.constants as scc
import  matplotlib.pyplot as plt

species = 'electrons'
fieldname = 'Ey'
lambda0 = .81e-6

# --- custom functions --- #
# ------------------------ #
omega0 = 2*np.pi*scc.c/lambda0
# Read field fieldname and return normalized max
def get_a0(res_dir, snapshot):
    header = res_dir + '/Header'
    print( snapshot )
    allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
    F = allrd[ fieldname ]
    return info['z'][-1], np.max(np.abs(F)) * scc.e/(scc.m_e*omega0*scc.c)
# Convert elements of a list to numpy arrays
def convert_to_np_array(list_in):
    list_out = []
    for elem in list_in:
        list_out.append( np.array( elem ) )
    return list_out

# --- MPI parallelization --- #
# --------------------------- #
# Get ordered list of snapshot files
res_dir = './lab_frame_data/';
file_list = glob.glob(res_dir + '/snapshot?????')
file_list.sort()
# Total number of files
nfiles = len(file_list)
# Each file has a unique number
number_list = range(nfiles)
comm_world = MPI.COMM_WORLD
me = comm_world.Get_rank()
nrank = comm_world.Get_size()
# List of files to process for current proc
my_list = file_list[ (me*nfiles)/nrank : ((me+1)*nfiles)/nrank ]
# List if file numbers for current proc
my_number_list = number_list[ (me*nfiles)/nrank : ((me+1)*nfiles)/nrank ]

# --- Run parallel analysis --- #
# ----------------------------- #
# Each MPI rank reads roughly (nb snapshots)/(nb ranks) snapshots.
# Works with any number of snapshots.
for count, filename in enumerate(my_list):
    zwin, a0 = get_a0( res_dir, filename )
    uzlab = read_raw_data.get_particle_field(filename, species, 'uz')/scc.c
    select_particles = (uzlab > 5.)
    uzlab = uzlab[ select_particles ]
    uzmean = np.mean(uzlab)

# --- gather and rank 0 plots --- #
# ------------------------------- #
# Gather particle quantities to rank 0 to plot history of quantities.
UZMEAN = comm_world.gather(uzmean, root=0)
ZWIN = comm_world.gather(zwin, root=0)
A0 = comm_world.gather(a0, root=0)
# Rank 0 does the plot.
if me == 0:
    # Convert to numpy arrays
    UZMEAN, ZWIN, A0 = convert_to_np_array([UZMEAN, ZWIN, A0])
    # Plot and save
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.plot(ZWIN*1.e3, UZMEAN)
    plt.xlabel('z (mm)')
    plt.ylabel('uz')
    plt.title( 'beam energy' )
    plt.grid()
    plt.subplot(2,1,2)
    plt.plot(ZWIN*1.e3, A0)
    plt.ylabel('a0')
    plt.title( 'beam propag angle' )
    plt.grid()
    fig.savefig('./image.png', bbox_inches='tight')
    plt.close()
