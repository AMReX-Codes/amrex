# Sphere of radius 0.75 centered at x = 1.0 in the fully periodic domain
# [-1.2,1.2]^3, together with its periodic image at x = -1.4, so that the body
# is consistent with the periodicity while crossing the x faces.  Ghost nodes
# sample the implicit function at their unwrapped coordinates, as in the
# legacy generator, which is only meaningful for a periodic body; the exported
# level set and cell flags must therefore agree on both sides of the seam.
# Fluid volume: 2.4^3 - 4/3*pi*0.75^3 = 12.0567.
nx = 32
ny = 32
nz = 32

max_grid_size = 16
algorithm_tests = 0

geometry.is_periodic = 1 1 1

eb2.geometry_method = marching_cubes
eb2.geom_type = parser
eb2.parser_function = "0.5625 - min((x-1.0)*(x-1.0) + y*y + z*z, (x+1.4)*(x+1.4) + y*y + z*z)"
