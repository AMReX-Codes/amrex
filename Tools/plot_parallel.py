import os
import glob
import matplotlib
import sys
import argparse
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
import matplotlib.pyplot as plt

'''
This script loops over all WarpX plotfiles in a directory and, for each
plotfile, saves an image showing the field and particles.

Requires yt>3.5 and Python3

It can be run serial:

> python plot_parallel.py --path <path/to/plt/files> --serial

or parallel

> mpirun -np 32 python plot_parallel.py --path <path/to/plt/files>

When running parallel, the plotfiles are distributed as evenly as possible
between MPI ranks.

This script also proposes an option to plot one quantity over all timesteps.
The data of all plotfiles are gathered to rank 0, and the quantity evolution
is plotted and saved to file. For the illustration, this quantity is the max
of Ey. Use " --plot_max_evolution Ey" to activate this option.

To get help, run
> python plot_parallel --help
'''

# Parse command line for options.
parser = argparse.ArgumentParser()
parser.add_argument('--path', default='diags/plotfiles',
                    help='path to plotfiles. Plotfiles names must be plt?????')
parser.add_argument('--plotlib', default='yt',
                    choices=['yt','matplotlib'],
                    help='Plotting library to use')
parser.add_argument('--field', default='Ez',
                    help='Which field to plot, e.g., Ez, By, jx or rho. The central slice in y is plotted')
parser.add_argument('--pjump', default=20,
                    help='When plotlib=matplotlib, we plot every pjump particle')
parser.add_argument('--vmax', type=float, default=None,
                    help='If specified, the colormap will have bounds [-vmax, vmax]')
parser.add_argument('--slicewidth', default=10.e-6,
                    help='Only particles with -slicewidth/2<y<slicewidth/2 are plotted')
parser.add_argument('--serial', action='store_true', default=False,
                    help='Specifies running in serial, avoiding the import of MPI')
parser.add_argument('--species', dest='pslist', nargs='+', type=str, default=None,
                    help='Species to be plotted, e.g., " --species beam plasma_e ". By default, all species in the simulation are shown')
parser.add_argument('--plot_max_evolution', type=str, default=None,
                    help='Quantity to plot the max of across all data files')
args = parser.parse_args()

plotlib = args.plotlib
plot_max_evolution = args.plot_max_evolution
vmax = args.vmax

# Sanity check
if int(sys.version[0]) != 3:
    print('WARNING: Parallel analysis was only tested with Python3')

matplotlib.rcParams.update({'font.size': 14})
pscolor = ['r','g','b','k','m','c','y','w']
pssize = 1.
yt_slicedir = {2:2, 3:1}
yt_aspect = {2:.05, 3:20}

# Get list of particle species.
def get_species(a_file_list):
    # if user-specified, just return the user list
    if args.pslist is not None:
        return args.pslist
    # otherwise, loop over all plotfiles to get particle species list
    psset = set()
    for filename in a_file_list:
        ds = yt.load( filename )
        for ps in ds.particle_types:
            if ps == 'all':
                continue
            psset.add(ps)
    pslist = list(psset)
    pslist.sort()
    return pslist

def plot_snapshot(filename):
    print( filename )
    # Load plotfile
    ds = yt.load( filename )
    # Get number of dimension
    dim = ds.dimensionality

    # Plot field colormap
    if plotlib == 'matplotlib':
        plt.figure(figsize=(12,7))
        # Read field quantities from yt dataset
        all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
        F = all_data_level_0['boxlib', args.field].v.squeeze()
        if dim == 3:
            F = F[:,int(F.shape[1]+.5)//2,:]
        extent = [ds.domain_left_edge[dim-1], ds.domain_right_edge[dim-1],
                  ds.domain_left_edge[0], ds.domain_right_edge[0]]
        # Plot field quantities with matplotlib
        plt.imshow(F, aspect='auto', extent=extent, origin='lower')
        plt.colorbar()
        plt.xlim(ds.domain_left_edge[dim-1], ds.domain_right_edge[dim-1])
        plt.ylim(ds.domain_left_edge[0], ds.domain_right_edge[0])
        if vmax is not None:
            plt.clim(-vmax, vmax)
    if plotlib == 'yt':
        # Directly plot with yt
        sl = yt.SlicePlot(ds, yt_slicedir[dim], args.field, aspect=yt_aspect[dim])
        if vmax is not None:
            sl.set_zlim(-vmax, vmax)

    # Plot particle quantities
    for ispecies, pspecies in enumerate(pslist):
        if pspecies in [x[0] for x in ds.field_list]:
            if plotlib == 'matplotlib':
                # Read particle quantities from yt dataset
                xp = all_data_level_0[pspecies, 'particle_position_x'].v
                if dim == 3:
                    yp = all_data_level_0[pspecies, 'particle_position_y'].v
                    zp = all_data_level_0[pspecies, 'particle_position_z'].v
                    select = yp**2<(args.slicewidth/2)**2
                    xp = xp[select] ; yp = yp[select] ; zp = zp[select]
                if dim == 2:
                    zp = all_data_level_0[pspecies, 'particle_position_y'].v
                # Select randomly one every pjump particles
                random_indices = np.random.choice(xp.shape[0], int(xp.shape[0]/args.pjump))
                if dim == 2:
                    xp=xp[random_indices] ; zp=zp[random_indices]
                if dim == 3:
                    xp=xp[random_indices] ; yp=yp[random_indices] ; zp=zp[random_indices]
                plt.scatter(zp,xp,c=pscolor[ispecies],s=pssize, linewidth=pssize,marker=',')
            if plotlib == 'yt':
                # Directly plot particles with yt
                sl.annotate_particles(width=(args.slicewidth, 'm'), p_size=pssize,
                                      ptype=pspecies, col=pscolor[ispecies])
    # Add labels to plot and save
    iteration = int(filename[-5:])
    if plotlib == 'matplotlib':
        plt.xlabel('z (m)')
        plt.ylabel('x (m)')
        plt.title('%s at iteration %d, time = %e s'%(args.field, iteration, ds.current_time))
        plt.savefig('plt_%s_%s_%05d.png'%(args.field, plotlib, iteration), bbox_inches='tight', dpi=300)
        plt.close()
    if plotlib == 'yt':
        sl.annotate_grids()
        sl.save('plt_%s_%s_%05d.png'%(args.field, plotlib, iteration))

# Compute max of field a_field in plotfile filename
def get_field_max( filename, a_field ):
    # Load plotfile
    ds = yt.load( filename )
    # Get number of dimension
    dim = ds.dimensionality
    # Read field quantities from yt dataset
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F = all_data_level_0['boxlib', a_field].v.squeeze()
    zwin = (ds.domain_left_edge[dim-1]+ds.domain_right_edge[dim-1])/2
    maxF = np.amax(F)
    return zwin, maxF

def plot_field_max(zwin_arr, maxF_arr):
    plt.figure()
    plt.plot(zwin_arr, maxF_arr)
    plt.xlabel('z (m)')
    plt.ylabel('%s (S.I.)'%plot_max_evolution)
    plt.title('Field max evolution')
    plt.savefig('max_%s_evolution.pdf'%plot_max_evolution, bbox_inches='tight')

### Analysis ###

# Get list of plotfiles
file_list = glob.glob(os.path.join(args.path, 'plt?????'))
file_list.sort()
nfiles = len(file_list)

# Get list of particle speciess to plot
pslist = get_species(file_list);

rank = 0
size = 1
if not args.serial:
    try:
        from mpi4py import MPI
        comm_world = MPI.COMM_WORLD
        rank = comm_world.Get_rank()
        size = comm_world.Get_size()
    except ImportError:
        pass

if rank == 0:
    print('number of MPI ranks: %d'%size)
    print('Number of plotfiles: %s'%nfiles)
    print('list of species: ', pslist)

if plot_max_evolution is not None:
    # Fill with a value less than any possible value
    zwin = np.full(nfiles, np.finfo(float).min)
    maxF = np.full(nfiles, np.finfo(float).min)

# Loop over files, splitting plotfile list among MPI ranks
# - plot field snapshot
# - store window position and field max in arrays
for count, filename in enumerate(file_list):
    if count%size != rank:
        continue

    plot_snapshot( filename )
    if plot_max_evolution is not None:
        zwin[count], maxF[count] = get_field_max( filename, plot_max_evolution )

if plot_max_evolution is not None:
    if size > 1:
        global_zwin = np.empty_like(zwin)
        global_maxF = np.empty_like(maxF)
        comm_world.Reduce(zwin, global_zwin, op=MPI.MAX)
        comm_world.Reduce(maxF, global_maxF, op=MPI.MAX)
        zwin = global_zwin
        maxF = global_maxF
    if rank == 0:
        plot_field_max(zwin, maxF)

