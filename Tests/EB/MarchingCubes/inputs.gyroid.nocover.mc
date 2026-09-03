# Coarsely resolved gyroid (parser implicit function).  Its saddles produce
# ambiguous MC33 cells whose fluid splits into more than one component, which
# the generator can only handle through the legacy-style nodal repair; with
# eb2.cover_multiple_cuts = 0 the build must abort with a message pointing at
# that parameter.
nx = 32
ny = 32
nz = 32

max_grid_size = 16
algorithm_tests = 0
custom_stl_test = 1

eb2.geometry_method = marching_cubes
eb2.geom_type = parser
eb2.parser_function = "sin(6*x)*cos(6*y)+sin(6*y)*cos(6*z)+sin(6*z)*cos(6*x)-0.9"
eb2.cover_multiple_cuts = 0
