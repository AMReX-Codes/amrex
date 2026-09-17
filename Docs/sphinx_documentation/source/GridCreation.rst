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
which is expensive for large domains.

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

Other applications do not want the grids decomposed in one coordinate
direction (for example, atmospheric codes that solve implicitly along vertical
columns).  Setting :cpp:`amr.no_chop_dir` to that direction (0 for *x*, 1 for
*y*, 2 for *z*; the default of -1 disables this) has the following effects.

- :cpp:`max_grid_size` and :cpp:`refine_grid_layout` are ignored in that
  direction on every level, and the :cpp:`blocking_factor` does not need to
  divide :cpp:`n_cell` in it.  (Applications built on :cpp:`class Amr` rather
  than :cpp:`AmrCore` still need an even :cpp:`n_cell` in every direction.)

- On level 0, if the :cpp:`blocking_factor` is 1 in the other directions, the
  domain is split in those directions only, into nearly equal grids that all
  span the domain in :cpp:`amr.no_chop_dir`.  Each direction gets the fewest
  pieces allowed by :cpp:`max_grid_size`.  If there are fewer grids than MPI
  processes and :cpp:`refine_grid_layout` permits, the number of pieces is
  doubled in the direction with the longest grids until there are enough.
  Any :cpp:`n_cell` works.  With a larger level 0 :cpp:`blocking_factor` the
  usual algorithm is used.

- On finer levels, the grids produced by the clustering are merged along that
  direction so that no two grids share an interior face normal to it, and
  they are never chopped in it afterwards.  A grid still covers only the part of the domain
  where cells are tagged; two tagged regions at different heights in the same
  column give two grids that do not touch.  Use
  :cpp:`amr.refine_whole_domain_dir` as well if every grid must span the entire
  domain in that direction.
  Grids at opposite ends of a periodic domain can still touch through the
  periodic boundary.

- The rules of :ref:`sec:grid_creation:odd` for domain sizes and
  :cpp:`max_grid_size` apply to every fine level, including levels with an
  even refinement ratio.

Users often like to ensure that coarse/fine boundaries are not too close to tagged cells; the
way to do this is to set :cpp:`amr.n_error_buf` to a large integer value (the default is 1).
This parameter is used to increase the number of tagged cells before the grids are defined;
if cell "*(i,j,k)*" satisfies the tagging criteria, then, for example, if :cpp:`amr.n_error_buf` is 3,
all cells in the 7x7x7 box from lower corner "*(i-3,j-3,k-3)*" to "*(i+3,j+3,k+3)*" will be tagged.

.. _sec:grid_creation:odd:

Odd Refinement Ratios and Odd Domain Sizes
------------------------------------------

Some applications, for example atmospheric codes that nest a fine domain
inside a coarse one, use odd refinement ratios such as 3 together with
domains whose sizes are not powers of 2, e.g., :cpp:`n_cell = 749 679 69`.
The usual advice of powers of 2 everywhere does not apply to them.  This
section explains how to choose :cpp:`blocking_factor` and
:cpp:`max_grid_size` in that situation.  It applies to the levels created
with an odd refinement ratio, and to all fine levels when
:cpp:`amr.no_chop_dir` is set.  Other levels follow the usual rules, even
when another level of the hierarchy has an odd refinement ratio.

**Blocking factor on the fine levels.**  On a level with refinement ratio
:math:`r`, choose a :cpp:`blocking_factor` that is :math:`r` times a power of
2, for example 24 (that is, :math:`3 \times 8`) for :math:`r = 3`, or a plain
power of 2 such as 8 in a direction where :math:`r = 1`.  The grids on that
level are then multiples of the :cpp:`blocking_factor`, and regridding
becomes much cheaper, because the grid generator works with blocks of that
many coarse cells in each direction instead of individual cells (with 8,
that is 512 times fewer cells in 3D).  A :cpp:`blocking_factor` of 1 is
allowed, but it makes regridding expensive on large domains and produces
many small grids, so it is best avoided.  A power of 2 that is not a
multiple of :math:`r`, such as 8 with :math:`r = 3`, is accepted for
backward compatibility, but the grids are then multiples of 6 rather than
8 (reported when :cpp:`amr.v` is positive).  In all cases
:cpp:`blocking_factor` divided by :math:`r`, rounded down, must be a power
of 2 (or less than 1), so 16 is rejected for :math:`r = 3`.

**Max grid size on the fine levels.**  :cpp:`max_grid_size` must be a
multiple of the :cpp:`blocking_factor`, e.g., 96 or 192 for a
:cpp:`blocking_factor` of 24, and at least twice the :cpp:`blocking_factor`
when :cpp:`n_cell` is not divisible by it.  Larger values give fewer, larger
grids; smaller values give more grids to distribute across processes.

**Level 0.**  The :cpp:`blocking_factor` on level 0 must divide
:cpp:`n_cell`, so when :cpp:`n_cell` has no convenient factors it has to be
1.  This costs nothing: the level 0 blocking factor has no effect on
regridding, and the level 0 grids do not need to line up with the blocking
factor of level 1.  :cpp:`max_grid_size` on level 0 can be anything; it just
sets how the domain is split among processes.

**Domain sizes.**  The domain does not need to be divisible by the
:cpp:`blocking_factor`.  In a non-periodic direction, the grids that reach
the upper domain boundary are simply cut off there, but no grid is ever
thinner than the :cpp:`blocking_factor`.  In a periodic direction, choose
:cpp:`n_cell` to be a multiple of the :cpp:`blocking_factor` divided by
:math:`r` (8 in the example above); if it is not, a smaller power of 2 is
used in that direction, the grids there are only multiples of :math:`r`
times that smaller number, and a warning is printed.

For example, with :cpp:`n_cell = 749 679 69`, :cpp:`ref_ratio_vect = 3 3 1`
and :cpp:`max_level = 1`, a good choice is

.. code-block:: none

   amr.blocking_factor_x = 1 24
   amr.blocking_factor_y = 1 24
   amr.blocking_factor_z = 1 8
   amr.max_grid_size_x   = 188 96
   amr.max_grid_size_y   = 188 96
   amr.max_grid_size_z   = 69 72

- The level 0 blocking factor is 1 because 749, 679 and 69 are odd, so no
  larger power of 2 divides them.
- On level 1 the blocking factor is :math:`3 \times 8 = 24` in *x* and *y*,
  where the refinement ratio is 3, and :math:`1 \times 8 = 8` in *z*, where it
  is 1.  Level 1 grids are therefore multiples of 24 by 24 by 8 cells, except
  where they are cut off at the upper domain boundary.
- The level 0 max grid size of 188 splits 749 cells into four grids in *x*
  and *y*; 69 in *z* keeps each level 0 grid whole in the vertical.
- The level 1 max grid size is 96, a multiple of 24, in *x* and *y*, and 72,
  a multiple of 8 that exceeds the 69 cells of the domain, in *z*, so that
  level 1 grids are never split in the vertical either.

Codes that require every grid to span whole vertical columns should also set
:cpp:`amr.no_chop_dir = 2`, which then takes care of the level 0 layout and of
the vertical extent of the fine grids automatically; see above.

