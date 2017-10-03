Input parameters
================

.. warning::

   This section is currently in development.


Overall simulation parameters
-----------------------------

* ``max_step`` (`integer`)
    The number of PIC cycles to perform.


Setting up the field mesh
-------------------------

* ``amr.n_cell`` (`2 integers in 2D`, `3 integers in 3D`)
    The number of grid points along each direction (on the **coarsest level**)

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

Distribution across MPI ranks and parallelization
-------------------------------------------------


* ``amr.max_grid_size`` (`integer`)
    Maximum allowable size of each **subdomain**
    (expressed in number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI rank.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI ranks
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.

* ``warpx.load_balance_int`` (`integer`)
    How often WarpX should try to redistribution the work across MPI ranks,
    in order to have better load balancing (expressed in number of PIC cycles
    inbetween two consecutive attempts at redistributing the work).
    Use 0 to disable load_balancing.

    When performing load balancing, WarpX measures the wall time for
    computational parts of the PIC cycle. It then uses this data to decide
    how to redistribute the subdomains across MPI ranks. (Each subdomain
    is unchanged, but its owner is changed in order to have better performance.)
    This relies on each MPI rank handling several (in fact many) subdomains
    (see ``max_grid_size``).


Particle initialization
-----------------------


Numerics and algorithms
-----------------------

* ``warpx.use_filter`` (`0 or 1`)
    Whether to smooth the charge and currents on the mesh, after depositing
    them from the macroparticles. This uses a bilinear filter
    (see the sub-section **Filtering** in :doc:`../theory/theory`).

* ``algo.current_deposition`` (`integer`)
    The algorithm to be used for current deposition:

     - ``0``: Esirkepov deposition, vectorized
     - ``1``: Esirkepov deposition, non-optimized
     - ``2``: Direct deposition, vectorized
     - ``3``: Direct deposition, non-optimized

* ``algo.charge_deposition`` (`integer`)
    The algorithm to be used for the charge density deposition:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

* ``algo.field_gathering`` (`integer`)
    The algorithm to be used for field gathering:

     - ``0``: Vectorized version
     - ``1``: Non-optimized version

* ``algo.particle_pusher`` (`integer`)
    The algorithm to be used for the particle pusher:

     - ``0``: Boris pusher
     - ``1``: Vay pusher

Diagnostics and output
----------------------

* ``amr.plot_int`` (`integer`)
    The number of PIC cycles inbetween two consecutive data dumps. Use a
    negative number to disable data dumping.
