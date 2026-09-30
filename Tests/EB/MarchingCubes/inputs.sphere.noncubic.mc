# Non-cubic cells (dz != dx): the marching-cubes generator must abort with a
# clear message instead of building an inconsistent geometry.
nx = 64
ny = 64
nz = 32

max_grid_size = 32
algorithm_tests = 0

eb2.geometry_method = marching_cubes
eb2.geom_type = sphere
eb2.sphere_center = 0.0 0.0 0.0
eb2.sphere_radius = 0.75
eb2.sphere_has_fluid_inside = 0
