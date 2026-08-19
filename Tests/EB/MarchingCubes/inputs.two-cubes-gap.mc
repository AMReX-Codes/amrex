nx = 16
ny = 16
nz = 16
max_grid_size = 8

algorithm_tests = 0
custom_stl_test = 1

geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 16.0 16.0 16.0

# Ordinary EB2::Build repairs the six multi-component cells instead of
# retaining multiple VoFs.
eb2.geometry_method = marching_cubes
eb2.stl_file = two_cubes_narrow_gap.stl
eb2.small_volfrac = 0.0
eb2.cover_multiple_cuts = 1
eb2.max_grid_size = 8
# eb2.mc_stl_file = two_cubes_narrow_gap_single_mc_triangles.stl
eb2.eb_surface_stl_file = two_cubes_narrow_gap_single_mc.stl
