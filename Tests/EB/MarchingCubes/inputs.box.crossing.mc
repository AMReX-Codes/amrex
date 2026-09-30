# Box that crosses the -y and +z domain faces with the default domain-face
# extension enabled: the geometry is extruded straight outward from the
# physical faces into the ghost region, and the in-domain volume is unchanged.
nx = 64
ny = 64
nz = 64

max_grid_size = 32
algorithm_tests = 0

eb2.geometry_method = marching_cubes
eb2.extend_domain_face = 1
eb2.geom_type = box
eb2.box_lo = -0.55 -1.5 -0.55
eb2.box_hi =  0.55  0.47  1.5
eb2.box_has_fluid_inside = 0
