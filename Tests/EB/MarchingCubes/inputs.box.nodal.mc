# Box faces lying exactly on grid nodes: exercises the exact-zero (ON) node
# handling of the marching-cubes generator with an implicit function.
nx = 64
ny = 64
nz = 64

xmin = -1.0
xmax = 1.0
ymin = -1.0
ymax = 1.0
zmin = -1.0
zmax = 1.0

max_grid_size = 32
algorithm_tests = 0

eb2.geometry_method = marching_cubes
eb2.geom_type = box
eb2.box_lo = -0.5 -0.5 -0.5
eb2.box_hi =  0.5  0.5  0.5
eb2.box_has_fluid_inside = 0
# The box faces lie on grid nodes, so cells owning a node-coincident EB patch
# legitimately have unit volume and a nonzero boundary area.
allow_full_volume_cut_cells = 1
