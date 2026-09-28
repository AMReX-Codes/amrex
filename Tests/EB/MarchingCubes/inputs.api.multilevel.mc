# EB2::Build(makeShop(SphereIF), Vector<Geometry>, ...) builds every coarse
# level directly with the marching-cubes generator instead of coarsening.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0
api_build = sphere
api_all_levels = 1
required_coarsening_level = 2
max_coarsening_level = 2

eb2.geometry_method = marching_cubes
eb2.sphere_center = 0.0 0.0 0.0
eb2.sphere_radius = 0.7
eb2.sphere_has_fluid_inside = 0
# The STL output must hold the finest level (same facet count as the
# single-level build in inputs.api.mc), not the last coarse level built.
eb2.mc_stl_file = api_multilevel_mc.stl
expected_mc_stl_facets = 13112
