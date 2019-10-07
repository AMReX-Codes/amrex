'''
This script loops over 3D plotfiles plt*****, generates a 3D rendering
of the data with fields and particles, and saves one image per plotfile to
plt_****_img.png. It was written for a laser-wakefield acceleration
simulation, and contains a lot of custom values (for transparency,
color intensity etc.), so feel free to modify it to meet your needs.

Execute the file with
> python video_yt.py
to generate the images. It can be quite slow for even moderately large
plotfiles.
'''

# Import statements
import sys, os
import yt, glob
import numpy as np

clight = 299792458.0

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
    box_patch = yt.visualization.volume_rendering.render_source.BoxSource(
        left_edge=ds.index.grids[1].LeftEdge+np.array([0., 0., z_shift])*yt.units.meter,
        right_edge=ds.index.grids[1].RightEdge+np.array([0., 0., z_shift])*yt.units.meter,
        color=[1.,0.1,0.1,.01])

# Handle 2 populations of particles: beam and plasma electrons
    # Color for each of these particles
    colors0_vect = [1., 1., 1., .05] # the last value is overwritten later
    colors1_vect = [1., 1., 1., .05] # the last value is overwritten later
    # particle0: read data and create a yt source object
    x0 = ad['particle0','particle_position_x'].v
    y0 = ad['particle0','particle_position_y'].v
    z0 = ad['particle0','particle_position_z'].v
    vertices0 = np.column_stack((x0,y0,z0))
    colors0 = np.tile(colors0_vect,(vertices0.shape[0], 1))
    colors0[:,3] = .01
    point0 = yt.visualization.volume_rendering.render_source.PointSource(vertices0, colors=colors0, radii=2)
    # particle1: read data and create a yt source object
    x1 = ad['particle1','particle_position_x'].v
    y1 = ad['particle1','particle_position_y'].v
    z1 = ad['particle1','particle_position_z'].v
    # select only some particles
    selector = np.abs(x1)<.1e-6
    x1 = x1[selector]
    y1 = y1[selector]
    z1 = z1[selector]
    vertices1 = np.column_stack((x1,y1,z1))
    colors1 = np.tile(colors1_vect,(vertices1.shape[0], 1))
    colors1[:,3] = .002
    point1 = yt.visualization.volume_rendering.render_source.PointSource(vertices1, colors=colors1, radii=1)

# Set the field rendering and camera attributes
    # Max field in the simulation. This is easy to get from a single image, but
    # it must be set by hand for a video
    my_max = 500000000000
    sc = yt.create_scene(ds, field='Ez')
    # Set camera properties
    cam = sc.camera
    cam.set_width(ds.quan(35, yt.units.micrometer))
    cam.position = ds.domain_center + np.array([15., 20., -5.])*yt.units.micrometer + np.array([0., 0., z_shift])*yt.units.meter
    cam.focus = ds.domain_center + np.array([0., 0., z_shift])*yt.units.meter
    cam.resolution = (2048, 2048)
    # Field rendering properties
    source = sc[0]
    source.set_field('Ez')
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
# Plot user-defined sources (here, 2 species for particles and 1 box)
    sc.add_source(point0)
    sc.add_source(point1)
    sc.add_source(box_patch)
# Save file
    sc.save(filename + '_img.png', sigma_clip=1.)
    return 0

# Get plt folders in current folder and loop over them.
file_list = glob.glob('./plt?????')
for filename in file_list:
    # disabled test: do not plot image if already exists
    # if os.path.isfile(filename + '.png') is False:
    print(filename)
    img_onestep(filename)
