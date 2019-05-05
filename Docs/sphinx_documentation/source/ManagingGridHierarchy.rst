.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _ss:grid_creation:

Grid Creation
-------------

To run an AMReX-based application you must specifiy the domain size by
specifying :cpp`n_cell` -- this is the number of cells spanning the domain 
in each coordinate direction at level 0.

Users often specify :cpp`max_grid_size` as well. The default load balancing algorithm then divides the 
domain in every direction so that each grid is no longer than :cpp`max_grid_size` in that direction.

Another popular input is :cpp`blocking_factor`.  The value of :cpp`blocking_factor` 
constrains grid creation in that in that each grid must be divisible by :cpp`blocking_factor`.  

Notes: 
 - :cpp`n_cell` must be given as three separate integers, one for each coordinate direction.
 - if :cpp`max_grid_size` is specified as a single integer *m* then each grid at each level must 
   be no larger than *m* in each coordinate direction
 - if :cpp`max_grid_size` is specified as multiple integers then the first 
   integer applies to level 0, the second to level 1, etc.  If you don't specify as many
   integers as there are levels, the final value will be used for the remaining levels.
 - if different values of :cpp`max_grid_size` are wanted for different coordinate directions, 
   then :cpp`max_grid_size_x`, :cpp`max_grid_size_y` and :cpp`max_grid_size_z` must be used.  
   If you don't specify as many integers as there are levels, the final value will be used for the remaining levels.
 - The value of :cpp`blocking_factor` also constrains the size of each grid in that all
   grids must be divisible by :cpp`blocking_factor`.  
 - The original domain (as specified by :cpp`n_cell`) must be divisible by :cpp`blocking_factor`
 - :cpp`max_grid_size` is not allowed to be less than :cpp`blocking_factor`.
 - to create identical grids of a specific size, e.g. of length *m* in each direction, 
   then set :cpp`max_grid_size` = *m* and :cpp`clocking_factor` = *m*.
 - note that :cpp`max_grid_size` is just an upper bound; if :cpp`n_cell = 48` 
   and :cpp`max_grid_size = 32`, then we will typically have one grid of length 32
   and one of length 16.

The gridding algorithm proceeds as follows:

#. If at level 0, the domain is initially defined by :cpp:`n_cell`
   as specified in the inputs file. If at level greater than 0,
   grids are created using the Berger-Rigoutsis clustering algorithm applied to the
   tagged cells from the section on :ref:`ss:regridding`, modified to ensure that
   all new fine grids are divisible by :cpp:`blocking_factor`.

#. Next, the grid list is chopped up if any grids are larger than :cpp:`max_grid_size`.
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

