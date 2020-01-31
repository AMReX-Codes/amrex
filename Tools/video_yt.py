#! /usr/bin/env python
# Copyright 2017-2020 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

'''
This script loops over 3D plotfiles plt*****, generates a 3D rendering
of the data with fields and particles, and saves one image per plotfile to
plt_****_img.png. It was written for a laser-wakefield acceleration
simulation, and contains a lot of custom values (for transparency,
color intensity etc.), so feel free to modify it to meet your needs.

Execute the file with
> python video_yt.py
or
> mpirun -np 4 python video_yt.py
to generate the images. It can be quite slow for even moderately large
plotfiles.
'''

# Import statements
import yt, glob
yt.enable_parallelism()
import numpy as np

field = 'Ez'
my_max = int(5.e9) # Field maximum amplitude
do_particles = True
species0 = 'beam'
species1 = 'electrons'
do_patch = False # if want to plot an MR patch
resolution = (512, 512)
camera_position = np.array([15., 20., -5.])*yt.units.micrometer
file_list = glob.glob('./diags/plotfiles/plt?????')

clight = 299792458.0 # must be the same value as in WarpX

def plot_species(species, ad, radii, transparency, abs_xmax):
    # Color for each of these particles
    colors_vect = [1., 1., 1., .05] # the last value is overwritten later
    x = ad[species,'particle_position_x'].v
    y = ad[species,'particle_position_y'].v
    z = ad[species,'particle_position_z'].v
    selector = np.abs(x) < abs_xmax
    x = x[selector] ; y = y[selector] ; z = z[selector]
    vertices = np.column_stack((x,y,z))
    colors = np.tile(colors_vect,(vertices.shape[0], 1))
    colors[:,3] = transparency
    point = yt.visualization.volume_rendering.render_source.PointSource(vertices, colors=colors, radii=radii)
    return point

# Create the 3d image for 1 timestep
# filename is the name of the folder (e.g. plt00000)
def img_onestep(filename):
    # Load the data
    ds = yt.load( filename )
    ad = ds.all_data()

    # Calculate the z position of the box.
    # You can use ds.domain_right_edge[2] instead. However, if a moving window
    # was used in the simulation, the rendering shows some jitter.
    # This is because a cell is added in z at some iterations but not all.
    # These lines calculate this jitter z_shift and remove it from the camera position and focus
    iteration=int(filename[-5:])
    dt = 1./clight * 1./np.sqrt((1./ad['dx'][-1]**2 + 1./ad['dy'][-1]**2 + 1./ad['dz'][-1]**2))
    z_front = dt * float(iteration) * clight
    z_shift = z_front-ds.domain_right_edge[2]

    # Create a yt source object for the level1 patch
    if do_patch:
        box_patch = yt.visualization.volume_rendering.render_source.BoxSource(
            left_edge=ds.index.grids[1].LeftEdge+np.array([0., 0., z_shift])*yt.units.meter,
            right_edge=ds.index.grids[1].RightEdge+np.array([0., 0., z_shift])*yt.units.meter,
            color=[1.,0.1,0.1,.01])

    # Handle 2 populations of particles: beam and plasma electrons
    if do_particles:
        point0 = plot_species(species0, ad, 2, .01, 1.)
        point1 = plot_species(species1, ad, 1, .002, 20.e-6)
    sc = yt.create_scene(ds, field=field)

    # Set camera properties
    cam = sc.camera
    dom_length = ds.domain_width[2].v
    cam.set_width(ds.quan(dom_length, yt.units.meter))
    cam.position = ds.domain_center + camera_position + np.array([0., 0., z_shift])*yt.units.meter
    cam.focus = ds.domain_center + np.array([0., 0., z_shift])*yt.units.meter
    cam.resolution = resolution
    # Field rendering properties
    source = sc[0]
    source.set_field(field)
    source.set_log(False)
    source.use_ghost_zones = True
    bounds = (-my_max, my_max)
    tf = yt.ColorTransferFunction(bounds)
    w = (.01*my_max)**2
    # Define the transfer function for 3d rendering
    # 3 isocontours for negative field values
    # The sharpness of the contour is controlled by argument width
    tf.add_gaussian(-.04 *my_max, width=8*w,  height=[0.1, 0.1, 1.0, 0.02])
    tf.add_gaussian(-.2 *my_max, width=5*w,  height=[0.1, 0.1, 1.0, 0.05])
    tf.add_gaussian(-.6 *my_max, width=w,  height=[0.0, 0.0, 1.0, 0.3])
    # 3 isocontours for positive field values
    tf.add_gaussian(.04 *my_max, width=8*w,  height=[1.0, 1.0, 0.2, 0.02])
    tf.add_gaussian(.2 *my_max, width=5*w,  height=[1.0, 1.0, 0.2, 0.05])
    tf.add_gaussian(.6 *my_max, width=w,  height=[1.0, 1.0, 0.0, 0.3])
    source.tfh.tf = tf
    source.tfh.bounds = bounds
    source.tfh.set_log(False)
    source.tfh.tf.grey_opacity = True
    if do_particles:
        sc.add_source(point0)
        sc.add_source(point1)
    if do_patch:
        sc.add_source(box_patch)
    sc.save('./img_' + filename[-8:] + '.png', sigma_clip=1.)

# Get plt folders in current folder and loop over them.
file_list.sort()
for filename in file_list[5:]:
    # disabled test: do not plot image if already exists
    # if os.path.isfile(filename + '.png') is False:
    print(filename)
    img_onestep(filename)
