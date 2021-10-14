amrex.fpe_trap_invalid = 1
eb_is_dirichlet = 0
eb_is_homog_dirichlet = 1

eb2.geom_type = channel

eb2.channel_askew = 1
eb2.channel_pt_on_top_wall = 2.0  0.0  0.0
eb2.channel_height = 1.0
eb2.channel_rotation = 10.0  60.0
eb2.channel_has_fluid_inside = 1

max_level = 0
ref_ratio = 2
n_cell = 80
max_grid_size =  40

max_coarsening_level = 0
prob_lo = 0.   0.   0.
prob_hi = 8.0  8.0  8.0
is_periodic = 0  0  0
scalars = 0  1

use_poiseuille = 1
poiseuille_askew = 1
poiseuille_pt_on_top_wall = 2.0  0.0  0.0
poiseuille_askew_rotation = 10.0  60.0
poiseuille_height = 1.0
poiseuille_no_flow_dir = 0

plot_file = "plot-3d-askew-all-mg"
