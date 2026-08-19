# Sphere crossing the +x face of a fully periodic domain.  The implicit
# function is not periodic, so the periodic images of the ghost nodes disagree
# with the function there; the generator must still build (root finding is
# bracketed with the function itself) instead of aborting.  No exact volume:
# the geometry itself is inconsistent with the periodicity.
nx = 32
ny = 32
nz = 32

max_grid_size = 16
algorithm_tests = 0
custom_stl_test = 1

geometry.is_periodic = 1 1 1

eb2.geometry_method = marching_cubes
eb2.geom_type = sphere
eb2.sphere_center = 1.0 0.2 0.3
eb2.sphere_radius = 0.75
eb2.sphere_has_fluid_inside = 0
