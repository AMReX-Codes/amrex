# Copyright 2019-2020 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import numpy as np

'''
This Python script helps a user to parallelize a WarpX simulation.

The user specifies the minimal size of the physical domain and the resolution
in each dimension, and the scripts computes:
- the number of cells and physical domain to satify the user-specified domain
  size and resolution AND make sure that the number of cells along each
  direction is a multiple of max_grid_size.
- a starting point on how to parallelize on Cori KNL (number of nodes, etc.).

When running in a boosted frame, the script also has the option to
automatically compute the number of cells in z to satisfy dx>dz in the boosted
frame.

Note that the script has no notion of blocking_factor. It is assumed that
blocking_factor = max_grid_size, and that all boxes have the same size.
'''

# Update the lines below for your simulation
# ------------------------------------------
# 2 elements for 2D, 3 elements for 3D
# Lower corner of the box
box_lo0 = np.array([-25.e-6, -25.e-6, -15.e-6])
# Upper corner of the box
box_hi0 = np.array([ 25.e-6,  25.e-6,  60.e-6])
# Cell size
dx = 1.e-6
dz = dx
cell_size = np.array([dx, dx, dz])
# Use this for simulations in a boosted frame if you
# want to enforce dz < dx / dx_over_dz_boosted_frame
compute_dz_boosted_frame = True
gamma_boost = 30.
dx_over_dz_boosted_frame = 1.1 # >1. is usually more stable
# ------------------------------------------

# similar to numpy.ceil, except the output data type is int
def intceil(num):
    return np.ceil(num).astype(int)

# Enlarge simulation boundaries to satisfy three conditions:
# - The resolution must be exactly the one provided by the user
# - The physical domain must cover the domain specified by box_lo0, box_hi0
# - The number of cells must be a multiple of mgs (max_grid_size).
def adjust_bounds(box_lo0, box_hi0, box_ncell0, mgs):
    cell_size = (box_hi0-box_lo0) / box_ncell0
    box_ncell = intceil(box_ncell0/mgs)*mgs
    box_lo = box_ncell * cell_size * box_lo0 / (box_hi0 - box_lo0)
    box_hi = box_ncell * cell_size * box_hi0 / (box_hi0 - box_lo0)
    return box_lo, box_hi, box_ncell

# Calculate parallelization for the simulation, given numerical parameters
# (number of cells, max_grid_size, number of threads per node etc.)
def nb_nodes_mpi(box_ncell,mgs,threadspernode,ompnumthreads,ngridpernode, ndim):
    nmpipernode = threadspernode/ompnumthreads
    ngridpermpi = ngridpernode/nmpipernode
    box_ngrids = box_ncell/mgs
    if ndim == 2:
        ngrids = box_ngrids[0] * box_ngrids[1]
    elif ndim == 3:
        ngrids = np.prod(box_ngrids)
    n_mpi = intceil( ngrids/ngridpermpi )
    n_node = intceil( n_mpi/nmpipernode )
    return n_node, n_mpi

# Get number of dimensions (2 or 3)
ndim = box_lo0.size
if compute_dz_boosted_frame:
    # Adjust dz so that dx/dz = dx_over_dz_boosted_frame in simulation frame
    cell_size[-1] = cell_size[0] / dx_over_dz_boosted_frame / 2. / gamma_boost
# Given the resolution, compute number of cells a priori
box_ncell0 = ( box_hi0 - box_lo0 ) / cell_size

if ndim == 2:
    # Set of parameters suitable for a 2D simulation on Cori KNL
    ngridpernode = 16.
    ompnumthreads = 8.
    mgs = 1024.
    threadspernode = 64. # HyperThreading level = 1: no hyperthreading
    distance_between_threads = int(68*4/threadspernode)
    c_option = int( ompnumthreads*distance_between_threads )
elif ndim == 3:
    # Set of parameters suitable for a 3D simulation on Cori KNL
    ngridpernode = 8.
    ompnumthreads = 8.
    mgs = 64.
    threadspernode = 64. # HyperThreading level = 1: no hyperthreading
    distance_between_threads = int(68*4/threadspernode)
    c_option = int( ompnumthreads*distance_between_threads )

# Adjust simulation bounds
box_lo, box_hi, box_ncell = adjust_bounds(box_lo0, box_hi0, box_ncell0, mgs)

# Calculate parallelization
n_node,n_mpi = nb_nodes_mpi(box_ncell, mgs, threadspernode, ompnumthreads, ngridpernode, ndim)

# Print results
string_output = ' ### Parameters used ### \n'
string_output += 'ngridpernode = ' + str(ngridpernode) + '\n'
string_output += 'ompnumthreads = ' + str(ompnumthreads) + '\n'
string_output += 'mgs (max_grid_size) = ' + str(mgs) + '\n'
string_output += 'threadspernode ( = # MPI ranks per node * OMP_NUM_THREADS) = ' + str(threadspernode) + '\n'
string_output += 'ndim = ' + str(ndim) + '\n\n'
string_output += 'box_lo = ' + str(box_lo) + '\n'
string_output += 'box_hi = ' + str(box_hi) + '\n'
string_output += 'box_ncell = ' + str(box_ncell) + '\n'
string_output += 'n_node = ' + str(n_node) + '\n'
string_output += 'n_mpi = ' + str(n_mpi) + '\n'
print(string_output)
