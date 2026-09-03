# Box that crosses the +x domain face, in a domain that is periodic in z, with
# the domain-face extension disabled.  Exercises periodic ghost handling
# and eb2.extend_domain_face = 0; only the part of the box
# inside the domain counts toward the expected volume.  (The body stays clear
# of the periodic faces so the geometry is consistent with the periodicity.)
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0

geometry.is_periodic = 0 0 1

eb2.geometry_method = marching_cubes
eb2.extend_domain_face = 0
eb2.geom_type = box
eb2.box_lo =  0.4 -0.7 -0.7
eb2.box_hi =  1.6  0.7  0.7
eb2.box_has_fluid_inside = 0
