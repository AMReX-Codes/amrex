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
eb2.stl_scale = 0.5

algorithm_tests = 0
# The STL faces lie on grid nodes, so cells owning a node-coincident EB patch
# legitimately have unit volume and a nonzero boundary area.
allow_full_volume_cut_cells = 1
