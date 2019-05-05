.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _ss:grid_creation:

Grid Creation
-------------

To run an AMReX-based application you must specifiy the domain size by
specifying :cpp:`n_cell` -- this is the number of cells spanning the domain 
in each coordinate direction at level 0.

Users often specify :cpp:`max_grid_size` as well. The default load balancing algorithm then divides the 
domain in every direction so that each grid is no longer than :cpp:`max_grid_size` in that direction.
If not specified by the user, :cpp:`max_grid_size` defaults to 128 in 2D and 32 in 3D (in each coordinate direction).

Another popular input is :cpp:`blocking_factor`.  The value of :cpp:`blocking_factor` 
constrains grid creation in that in that each grid must be divisible by :cpp:`blocking_factor`.  
Note that both the domain (at each level) and :cpp:`max_grid_size` must be divisible by :cpp:`blocking_factor`
If not specified by the user, :cpp:`blocking_factor` defaults to 8 in each coordinate direction.

There is one more default behavior to be aware of.  There is a boolean :cpp:`refine_grid_layout` 
that defaults to true but can be over-ridden at run-time.
If :cpp:`refine_grid_layout` is true and the number of grids created is less than the number of processors 
(Ngrids < Nprocs), then grids will be further subdivided until Ngrids >= Nprocs.

Caveat: if subdividing the grids to achieve Ngrids >= Nprocs would violate the 
:cpp:`blocking_factor` criterion then additional grids are not created and the 
number of grids will remain less than the number of processors

Note that :cpp:`n_cell` must be given as three separate integers, one for each coordinate direction.

However, :cpp:`max_grid_size` and :cpp:`blocking_factor` can be specified as a single value 
applying to all coordinate directions, or as separate values for each direction.  

 - if :cpp:`max_grid_size` (or :cpp:`blocking_factor`) is specified as multiple integers then the first 
   integer applies to level 0, the second to level 1, etc.  If you don't specify as many
   integers as there are levels, the final value will be used for the remaining levels.

 - if different values of :cpp:`max_grid_size` (or :cpp:`blocking_factor`) are wanted for each coordinate direction, 
   then :cpp:`max_grid_size_x`, :cpp:`max_grid_size_y` and :cpp:`max_grid_size_z` 
   (or :cpp:`blocking_factor_x`, :cpp:`blocking_factor_y` and :cpp:`blocking_factor_z`) must be used.  
   If you don't specify as many integers as there are levels, the final value will be used for the remaining levels.

Additional notes:

 - to create identical grids of a specific size, e.g. of length *m* in each direction, 
   then set :cpp:`max_grid_size` = *m* and :cpp:`blocking_factor` = *m*.

 - note that :cpp:`max_grid_size` is just an upper bound; with :cpp:`n_cell = 48` 
   and :cpp:`max_grid_size = 32`, we will typically have one grid of length 32 and one of length 16.

The gridding algorithm proceeds as follows:

#. If at level 0, the domain is initially defined by :cpp:`n_cell`
   as specified in the inputs file. If at level greater than 0,
   grids are created using the Berger-Rigoutsis clustering algorithm applied to the
   tagged cells from the section on :ref:`ss:regridding`, modified to ensure that
   all new fine grids are divisible by :cpp:`blocking_factor`.

#. Next, the grid list is chopped up if any grids are larger than :cpp::`max_grid_size`.
   Note that because :cpp:`max_grid_size` is a multiple of :cpp:`blocking_factor`
   (as long as :cpp:`max_grid_size` is greater than :cpp:`blocking_factor`),
   the blocking_factor criterion is still satisfied.

#. Next, if ``refine_grid_layout = 1`` and there are more processors than grids
   at this level, then the grids at this level are further divided in order to ensure that
   no processors has less than one grid (at each level).
   In :cpp:`AmrMesh::ChopGrids`,

   -  if :cpp:`max_grid_size / 2` in the :cpp:`BL_SPACEDIM` direction is a multiple of
      :cpp:`blocking_factor`, then chop the grids in the :cpp:`BL_SPACEDIM` direction
      so that none of the grids are longer in that direction than :cpp:`max_grid_size / 2`

   -  If there are still fewer grids than processes, repeat the procedure in the
      :cpp:`BL_SPACEDIM-1` direction, and again in the :cpp:`BL_SPACEDIM-2` direction if necessary

   -  If after completing a sweep in all coordinate directions with :cpp:`max_grid_size / 2`,
      there are still fewer grids than processes, repeat the steps above with :cpp:`max_grid_size / 4`.

