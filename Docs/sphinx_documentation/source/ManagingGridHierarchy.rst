.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _ss:managinggrid:

Managing the Grid Hierarchy
===========================

.. _ss:grid_creation:

Grid Creation
-------------

The gridding algorithm proceeds in this order, using the parameters described
in the section on the :ref:`ss:amrcore`.

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

