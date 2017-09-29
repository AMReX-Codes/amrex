Input parameters
================

Overall simulation parameters
-----------------------------

* ``max_step`` (`integer`)
    The number of PIC cycles to perform.


Setting up the field mesh
-------------------------

* ``amr.n_cell`` (`2 integers in 2D`, `3 integers in 3D`)
    The number of grid points along each direction (on the **coarsest level**)

* ``amr.max_grid_size`` (`integer`)
    Maximum allowable size of each subdomain
    (in terms of number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI process.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI processes
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.

* ``amr.max_level`` (`integer`)
    When using mesh refinement, the number of refinement levels that will be used.

    Use 0 in order to disable mesh refinement.

* ``geometry.is_periodic`` (`2 integers in 2D`, `3 integers in 3D`)
    Whether the boundary conditions are periodic, in each direction.

    For each direction, use 1 for periodic conditions, 0 otherwise.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters)
    The extent of the full simulation box. This box is rectangular, and thus its
    extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and
    upper corner (``geometry.prob_hi``).

* ``warpx.fine_tag_lo`` and ``warpx.fine_tag_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters)
    **When using static mesh refinement with 1 level**, the extent of the refined patch.
    This patch is rectangular, and thus its extent is given here by the coordinates
    of the lower corner (``warpx.fine_tag_lo``) and upper corner (``warpx.fine_tag_hi``).


Particle initialization
-----------------------


Numerics and algorithms
-----------------------

* ``warpx.use_filter`` (`0 or 1`)
    Whether to smooth the charge and currents on the mesh, after depositing
    them from the macroparticles. This uses a bilinear filter
    (see the sub-section **Filtering** in :doc:`../theory/theory`).

Diagnostics and output
----------------------

* ``amr.plot_int`` (`integer`)
    The number of PIC cycles inbetween two consecutive data dumps.
