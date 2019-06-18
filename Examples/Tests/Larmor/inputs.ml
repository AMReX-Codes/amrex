# Maximum number of time steps
max_step = 400

# number of grid points
amr.n_cell =  32  32

# The lo and hi ends of grids are multipliers of blocking factor
amr.blocking_factor = 16

# Maximum allowable size of each subdomain in the problem domain;
#    this is used to decompose the domain for parallel calculations.
amr.max_grid_size = 64

# Maximum level in hierarchy (for now must be 0, i.e., one level in total)
amr.max_level = 1

warpx.fine_tag_lo = -0.8  -0.8
warpx.fine_tag_hi =  0.8   0.8

amr.plot_int = 2   # How often to write plotfiles.  "<= 0" means no plotfiles.

warpx.plot_raw_fields = 1
warpx.plot_dive = 1
warpx.plot_divb = 1
warpx.plot_finepatch = 1
warpx.plot_crsepatch = 1

# Geometry
geometry.coord_sys   = 0                  # 0: Cartesian
geometry.is_periodic = 0     0     0      # Is periodic?
geometry.prob_lo     = -2.0  -2.0         # physical domain
geometry.prob_hi     =  2.0   2.0

# PML
warpx.do_pml = 1
warpx.pml_ncell = 10

warpx.B_external = 0.0  0.00078110417851950768  0.0

# Verbosity
warpx.verbose = 1

# Algorithms

# CFL
warpx.cfl = 1.0

# particles
particles.nspecies = 2
particles.species_names = electron positron

particles.nspecies = 1
particles.species_names = electron

electron.charge = -q_e
electron.mass = m_e
electron.injection_style = "SingleParticle"
electron.single_particle_pos = 0.0  0.0  -1.25
electron.single_particle_vel = -0.45825756949558416  0.0  0.0   # gamma*beta

positron.charge = q_e
positron.mass = m_e
positron.injection_style = "SingleParticle"
positron.single_particle_pos = 0.0  0.0  -1.25
positron.single_particle_vel = 0.45825756949558416  0.0  0.0   # gamma*beta

electron.single_particle_weight = 1.0e12
positron.single_particle_weight = 1.0e12

# interpolation
interpolation.nox = 3
interpolation.noy = 3
interpolation.noz = 3

# Moving window
warpx.do_moving_window = 0

warpx.do_dive_cleaning = 1
