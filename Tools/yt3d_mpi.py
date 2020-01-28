# Copyright 2018-2019 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

'''
This script loops over 3D plotfiles plt*****, generates a 3D rendering
of the data with fields and particles, and saves one image per plotfile to
img_*****.png. It was written for a beam-driven wakefield acceleration
simulation, and contains a lot of custom values (for transparency,
color intensity etc.), so feel free to modify it to meet your needs.

Execute the file with e.g.
> mpirun -np 12 python yt3d_mpi.py
to generate the images. It can be quite slow for even moderately large
plotfiles.
'''

import yt, glob
from mpi4py import MPI
import numpy as np
import scipy.constants as scc
yt.funcs.mylog.setLevel(50)

# my_max = 1.e11 # for smooth rendering
my_max = 5.e10 # for layered rendering
species_to_plot = ['plasma_e', 'beam', 'driver']
# For each species, provide [red, green, blue, alpha] between 0. and 1.
species_colors = { 'plasma_e': [1., 1., 1., .15],
                   'beam'    : [1., 1., 1., .2 ],
                   'driver'  : [1., 1., 1., .2 ] }
# provide these to avoid jitter when using a moving window
use_moving_window = True
plot_mr_patch = False
rendering_type = 'layers' # 'layers' or 'smooth'
maxwell_solver = 'ckc' # 'ckc' or 'yee'
cfl = 0.99
file_list = glob.glob('plotfiles/plt?????')

bounds = ( -my_max, my_max )
z_shift = 0.
w = (.01*my_max)**2

def jitter_shift(ds, ad, cfl, iteration):
    if maxwell_solver == 'yee':
        dt = 1./scc.c * 1./np.sqrt((1./ad['dx'][-1]**2 + 1./ad['dy'][-1]**2 + 1./ad['dz'][-1]**2))
    elif maxwell_solver == 'ckc':
        dt = cfl * min( [ ad['dx'][-1], ad['dy'][-1], ad['dz'][-1] ] )  / scc.c
    z_front = dt * float(iteration) * scc.c + 7.5e-6*yt.units.meter
    z_shift = z_front-ds.domain_right_edge[2]
    return z_shift

def get_species_ytpoints(ad, species, color_vec):
    xp = ad[species,'particle_position_x'].v
    yp = ad[species,'particle_position_y'].v
    zp = ad[species,'particle_position_z'].v
    if species == 'plasma_e':
        selection = np.abs(xp)<2.e-6
        zp = zp[selection]
        yp = yp[selection]
        xp = xp[selection]
    vertices = np.column_stack((xp,yp,zp))
    colors = np.tile(color_vec,(vertices.shape[0], 1))
    points = yt.visualization.volume_rendering.render_source.PointSource(vertices, colors=colors, radii=1)
    return points

def img_onestep(filename):
    ds = yt.load( filename )
    ad = ds.all_data()
    iteration=int(filename[-5:])
    sc = yt.create_scene(ds, field='Ez')
    if use_moving_window:
        z_shift = jitter_shift( ds, ad, cfl, iteration )
    array_shift = z_shift * np.array([0., 0., 1.])
    if plot_mr_patch:
        box_patch = yt.visualization.volume_rendering.render_source.BoxSource(
            left_edge =ds.index.grids[1].LeftEdge +array_shift,
            right_edge=ds.index.grids[1].RightEdge+array_shift,
            color=[1.,0.1,0.1,.01] )
        sc.add_source(box_patch)
    ########################
    ### volume rendering ###
    ########################
    source = sc[0]
    source.use_ghost_zones = True
    source.grey_opacity = True
    source.set_log(False)
    tf = yt.ColorTransferFunction(bounds)
    if rendering_type == 'smooth':
        tf.add_gaussian(-my_max/4, width=15**2*w,  height=[0.0, 0.0, 1.0, 1])
        tf.add_gaussian( my_max/4, width=15**2*w,  height=[1.0, 0.0, 0.0, 1])
    if rendering_type == 'layers':
        # NEGATIVE
        tf.add_gaussian(-.04 *my_max, width=8*w,  height=[0.1, 0.1, 1.0, 0.2])
        tf.add_gaussian(-.2 *my_max, width=5*w,  height=[0.1, 0.1, 1.0, 0.5])
        tf.add_gaussian(-.6 *my_max, width=w,  height=[0.0, 0.0, 1.0, 1.])
        # POSITIVE
        tf.add_gaussian(.04 *my_max, width=8*w,  height=[1.0, 1.0, 0.2, 0.2])
        tf.add_gaussian(.2 *my_max, width=5*w,  height=[1.0, 1.0, 0.2, 0.5])
        tf.add_gaussian(.6 *my_max, width=w,  height=[1.0, 1.0, 0.0, 1.])

    ######################
    ### plot particles ###
    ######################
    species_points = {}
    for species in species_to_plot:
        species_points[ species ] = get_species_ytpoints(ad,
                                        species, species_colors[species])
        sc.add_source( species_points[ species ] )
    source.tfh.tf = tf
    source.tfh.bounds = bounds
    #########################
    ### camera properties ###
    #########################
    cam = sc.camera
    cam.resolution = (2048, 2048)
    cam.width = .00018*yt.units.meter
    cam.focus = ds.domain_center + \
                np.array([0., 0., 10.e-6 ])*yt.units.meter + \
                array_shift
    cam.position = ds.domain_center + \
                np.array([15., 15., -5.  ])*yt.units.micrometer + \
                array_shift
    cam.normal_vector = [-0.3, -0.3, -.2]
    cam.switch_orientation()
    # save image
    if rendering_type == 'smooth':
        sc.save('img_' + str(my_number_list[count]).zfill(5), sigma_clip=5.)
    if rendering_type == 'layers':
        sc.save('img_' + str(my_number_list[count]).zfill(5), sigma_clip=2.)

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
for count, filename in enumerate(my_list):
    print('processing ' + filename)
    img_onestep(filename)
