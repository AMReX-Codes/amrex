# Coarse levels rebuilt from the implicit function (not coarsened) so that the
# marching-cubes generator runs through the GeometryShop rebuild path.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0

required_coarsening_level = 2
max_coarsening_level = 2
build_coarse_level_by_coarsening = 0

eb2.geometry_method = marching_cubes
eb2.geom_type = sphere
eb2.sphere_center = 0.0 0.0 0.0
eb2.sphere_radius = 0.7
eb2.sphere_has_fluid_inside = 0
