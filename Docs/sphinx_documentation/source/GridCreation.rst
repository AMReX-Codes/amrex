.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:grid_creation:

Grid Creation
-------------

To run an AMReX-based application you must specify the domain size by
specifying :cpp:`n_cell` -- this is the number of cells spanning the domain
in each coordinate direction at level 0.

Users often specify :cpp:`max_grid_size` as well. The default load balancing algorithm then divides the
domain in every direction so that each grid is no longer than :cpp:`max_grid_size` in that direction.
If not specified by the user, :cpp:`max_grid_size` defaults to 128 in 2D and 32 in 3D (in each coordinate direction).

Another popular input is :cpp:`blocking_factor`.  The value of :cpp:`blocking_factor`
constrains grid creation in that each grid must be divisible by :cpp:`blocking_factor`.
Note that both the domain (at each level) and :cpp:`max_grid_size` must be divisible by :cpp:`blocking_factor`,
and that :cpp:`blocking_factor` must be either 1 or a power of 2 (otherwise the gridding algorithm
would not in fact create grids divisible by  :cpp:`blocking_factor` because of how  :cpp:`blocking_factor`
is used in the gridding algorithm).  See :ref:`sec:grid_creation:odd` for the
exceptions that apply with odd refinement ratios.

If not specified by the user, :cpp:`blocking_factor` defaults to 8 in each coordinate direction.
The typical purpose of :cpp:`blocking_factor` is to ensure that the grids will be
sufficiently coarsenable for good multigrid performance.

The :cpp:`blocking_factor` on a level :math:`\ell > 0` also controls the cost of
regridding.  Before the tagged cells on level :math:`\ell-1` are clustered into
grids, they are coarsened by :cpp:`blocking_factor` on level :math:`\ell` divided
by the refinement ratio between the two levels.  A :cpp:`blocking_factor` of 1
therefore means that the clustering algorithm works on every tagged cell,
which is expensive for large domains.  The :cpp:`blocking_factor` on level 0 has
no effect on regridding; it only constrains how the level 0 grids are subdivided.

There is one more default behavior to be aware of.  There is a boolean :cpp:`refine_grid_layout`
that defaults to true but can be overridden at run time.
If :cpp:`refine_grid_layout` is true and the number of grids created is less than the number of processors
(Ngrids < Nprocs), then grids will be further subdivided until Ngrids >= Nprocs.

Caveat: if subdividing the grids to achieve Ngrids >= Nprocs would violate the
:cpp:`blocking_factor` criterion, then additional grids are not created and the
number of grids will remain less than the number of processors.

Note that :cpp:`n_cell` must be given as three separate integers, one for each coordinate direction.

However, :cpp:`max_grid_size` and :cpp:`blocking_factor` can be specified as a single value
applying to all coordinate directions, or as separate values for each direction.

 - If :cpp:`max_grid_size` (or :cpp:`blocking_factor`) is specified as multiple integers then the first
   integer applies to level 0, the second to level 1, etc.  If you don't specify as many
   integers as there are levels, the final value will be used for the remaining levels.

 - If different values of :cpp:`max_grid_size` (or :cpp:`blocking_factor`) are wanted for each coordinate direction,
   then :cpp:`max_grid_size_x`, :cpp:`max_grid_size_y` and :cpp:`max_grid_size_z`
   (or :cpp:`blocking_factor_x`, :cpp:`blocking_factor_y` and :cpp:`blocking_factor_z`) must be used.
   If you don't specify as many integers as there are levels, the final value will be used for the remaining levels.

Additional notes:

 - To create identical grids of a specific size, e.g. of length *m* in each direction,
   then set :cpp:`max_grid_size` = *m* and :cpp:`blocking_factor` = *m*.

 - Note that :cpp:`max_grid_size` is just an upper bound; with :cpp:`n_cell = 48`
   and :cpp:`max_grid_size = 32`, we will typically have one grid of length 32 and one of length 16.

The grid creation process at level 0 proceeds as follows (if not using the KD-tree approach):

#. The domain is initially defined by a single grid of size :cpp:`n_cell`.

#. If :cpp:`n_cell` is greater than :cpp:`max_grid_size` then the grids are subdivided until
   each grid is no longer than  :cpp:`max_grid_size` cells on each side.  The :cpp:`blocking_factor` criterion
   (i.e., that the length of each side of each grid is divisible by :cpp:`blocking_factor` in that direction)
   is satisfied during this process.

#. Next, if :cpp:`refine_grid_layout = true` and there are more processors than grids
   at this level, then the grids at this level are further divided until Ngrids >= Nprocs
   (unless doing so would violate the :cpp:`blocking_factor` criterion).

The creation of grids at levels > 0 begins by tagging cells at the coarser level and follows
the Berger-Rigoutsos clustering algorithm with the additional constraints of satisfying
the :cpp:`blocking_factor` and :cpp:`max_grid_size` criteria.  An additional parameter
becomes relevant here: the "grid efficiency", specified as :cpp:`amr.grid_eff` in the inputs file.
This threshold value, which defaults to 0.7 (or 70%), is used to ensure that
grids do not contain too large a fraction of un-tagged cells.   We note that the grid creation
process attempts to satisfy the :cpp:`amr.grid_eff` constraint but will not do so if it means
violating the :cpp:`blocking_factor` criterion.

Some applications want the fine levels to cover the entire domain in one coordinate
direction, no matter where the cells are tagged.  Setting :cpp:`amr.refine_whole_domain_dir`
to that direction (0 for *x*, 1 for *y*, 2 for *z*; the default of -1 disables this)
makes the tagging of a cell behave as if the whole line of cells through it in that
direction were tagged.  The clustering is then performed in one fewer dimension, so
:cpp:`amr.grid_eff` refers to the fraction of tagged cells in the plane perpendicular
to that direction.  The resulting grids may still be chopped in that direction by
:cpp:`max_grid_size` and :cpp:`refine_grid_layout`, but together they always cover the
entire domain.

Users often like to ensure that coarse/fine boundaries are not too close to tagged cells; the
way to do this is to set :cpp:`amr.n_error_buf` to a large integer value (the default is 1).
This parameter is used to increase the number of tagged cells before the grids are defined;
if cell "*(i,j,k)*" satisfies the tagging criteria, then, for example, if :cpp:`amr.n_error_buf` is 3,
all cells in the 7x7x7 box from lower corner "*(i-3,j-3,k-3)*" to "*(i+3,j+3,k+3)*" will be tagged.

.. _sec:grid_creation:odd:

Odd Refinement Ratios and Odd Domain Sizes
------------------------------------------

Some applications (for example, atmospheric codes that nest a fine domain
inside a coarse one) use odd refinement ratios such as 3 together with
domains whose sizes are not powers of 2, e.g., :cpp:`n_cell = 749 679 69`.
The rules above are relaxed for them as follows.  Consider a fine level
:math:`\ell > 0` with refinement ratio :math:`r` (which may differ by
direction) between levels :math:`\ell-1` and :math:`\ell`.

- The :cpp:`blocking_factor` on level :math:`\ell` may be :math:`r` times a
  power of 2 (e.g., 24 for :math:`r = 3`).  The tags on level :math:`\ell-1`
  are then coarsened by that power of 2 before clustering, and the grids on
  level :math:`\ell` are multiples of the :cpp:`blocking_factor` and are
  coarsenable by :math:`r`.  As always, :cpp:`max_grid_size` on level
  :math:`\ell` must be a multiple of the :cpp:`blocking_factor`.

- If the :cpp:`blocking_factor` on level :math:`\ell` is a power of 2 that is
  not a multiple of :math:`r` (e.g., 8 with :math:`r = 3`), the grids are
  multiples of :math:`\max(1, b/r) \cdot r` instead, where :math:`b` is the
  :cpp:`blocking_factor` (6 in the example).  A warning is printed.

- The :cpp:`blocking_factor` on level 0 must still divide :cpp:`n_cell`, so it
  is often 1 when :cpp:`n_cell` has no small factors.  This is harmless: it does
  not affect regridding, only how :cpp:`refine_grid_layout` may subdivide the
  level 0 grids.  The level 0 grids do not need to be aligned with the
  blocking factor of level 1.

- The domain does not need to be divisible by the level :math:`\ell`
  :cpp:`blocking_factor` in non-periodic directions.  The grids touching the
  upper domain boundary are simply truncated there.  In periodic directions,
  the level :math:`\ell-1` domain must be divisible by the
  :cpp:`blocking_factor` on level :math:`\ell` divided by :math:`r`.

For example, with :cpp:`n_cell = 749 679 69`, :cpp:`ref_ratio_vect = 3 3 1` and
:cpp:`max_level = 1`, the following gives level 1 grids that are multiples of
24 in *x* and *y* and coarsens the tags by 8 in every direction before
clustering:

.. code-block:: none

   amr.blocking_factor_x = 1 24
   amr.blocking_factor_y = 1 24
   amr.blocking_factor_z = 1 8
   amr.max_grid_size_x   = 188 96
   amr.max_grid_size_y   = 188 96
   amr.max_grid_size_z   = 69 72

