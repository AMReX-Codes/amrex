# Parser implicit function through the marching-cubes generator: a paraboloid.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0

eb2.geometry_method = marching_cubes
eb2.geom_type = parser
eb2.parser_function = "0.36 - (x*x + y*y) - 0.5*z"
