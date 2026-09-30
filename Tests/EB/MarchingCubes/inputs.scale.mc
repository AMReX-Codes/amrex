# Scale and decomposition invariance of the marching-cubes generator: the same two-sphere
# implicit function is built with the function multiplied by 1e-6, 1 and 1e6
# and with two different eb2.max_grid_size values, and must give the same fluid
# volume and nodal repair.  The main build is a
# plain sphere.
nx = 32
ny = 32
nz = 32

max_grid_size = 16
algorithm_tests = 0
scale_invariance_test = 1

eb2.geometry_method = marching_cubes
eb2.geom_type = sphere
eb2.sphere_center = 0.0 0.0 0.0
eb2.sphere_radius = 0.75
eb2.sphere_has_fluid_inside = 0
