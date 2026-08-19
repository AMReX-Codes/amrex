# EB2::Build(makeShop(SphereIF), ...) through the C++ API with the
# marching-cubes generator selected by eb2.geometry_method.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0
api_build = sphere

eb2.geometry_method = marching_cubes
eb2.sphere_center = 0.0 0.0 0.0
eb2.sphere_radius = 0.7
eb2.sphere_has_fluid_inside = 0
