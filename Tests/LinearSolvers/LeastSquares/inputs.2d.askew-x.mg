amrex.fpe_trap_invalid = 1
eb_is_dirichlet = 0
eb_is_homog_dirichlet = 1

eb2.geom_type = channel

eb2.channel_pt_on_top_wall = 0.0  1.8  0.0
eb2.channel_height = 1.0
eb2.channel_rotation = -30.
eb2.channel_has_fluid_inside = 1

max_level = 0
ref_ratio = 2
n_cell = 20
max_grid_size = 10
max_coarsening_level = 0
prob_lo =  0.   0.   0.
prob_hi = 2.0  2.0   0.
is_periodic = 0  0  0
scalars = 0  1

use_poiseuille = 1
poiseuille_pt_on_top_wall = 0.0  1.8  0.0
poiseuille_height = 1.0
poiseuille_rotation = -30.

plot_file = "plot-askew-x-mg"
