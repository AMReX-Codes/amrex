# Analytic box implicit function through the marching-cubes generator.
# Same body as cube.stl, so the expected volume matches inputs.cube.mc.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0

eb2.geometry_method = marching_cubes
eb2.geom_type = box
eb2.box_lo = -1.0 -1.0 -1.0
eb2.box_hi =  1.0  1.0  1.0
eb2.box_has_fluid_inside = 0
