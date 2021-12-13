.. _Chap:InputsLoadBalancing:

Gridding and Load Balancing
===========================

The following inputs must be preceded by "amr" and determine how we create the grids and how often we regrid.

+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| Parameter              | Description                                                           |   Type      | Default   |
+========================+=======================================================================+=============+===========+
| regrid_int             | How often to regrid (in number of steps at level 0)                   |   Int       |    -1     |
|                        | if regrid_int = -1 then no regridding will occur                      |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_x        | Maximum number of cells at level 0 in each grid in x-direction        |    Int      | 32        |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_y        | Maximum number of cells at level 0 in each grid in y-direction        |    Int      | 32        |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_z        | Maximum number of cells at level 0 in each grid in z-direction        |    Int      | 32        |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_x      | Each grid must be divisible by blocking_factor_x in x-direction       |    Int      |  8        |
|                        | (must be 1 or power of 2)                                             |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_y      | Each grid must be divisible by blocking_factor_y in y-direction       |    Int      |  8        |
|                        | (must be 1 or power of 2)                                             |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| blocking_factor_z      | Each grid must be divisible by blocking_factor_z in z-direction       |    Int      |  8        |
|                        | (must be 1 or power of 2)                                             |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| refine_grid_layout     | Split grids in half until the number of grids is no less than the     |    Bool     |  true     |
|                        | number of procs. (Will be overridden if refine_grid_layout_[x,y,z]    |             |           |
|                        | is specified)                                                         |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| refine_grid_layout_x   | Allow grids to be split in the x-dimension when refining the layout.  |    Int      |  1        |
|                        | (1 to allow or 0 to disallow)                                         |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| refine_grid_layout_y   | Allow grids to be split in the y-dimension when refining the layout.  |    Int      |  1        |
|                        | (1 to allow or 0 to disallow)                                         |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+
| refine_grid_layout_z   | Allow grids to be split in the z-dimension when refining the layout.  |    Int      |  1        |
|                        | (1 to allow or 0 to disallow)                                         |             |           |
+------------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "particles".

+-------------------+-----------------------------------------------------------------------+-------------+-----------+
|  Parameter        | Description                                                           |   Type      | Default   |
+===================+=======================================================================+=============+===========+
| max_grid_size_x   | Maximum number of cells at level 0 in each grid in x-direction        |    Int      | 32        |
|                   | for grids in the ParticleBoxArray if dual_grid is true                |             |           |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_y   | Maximum number of cells at level 0 in each grid in y-direction        |    Int      | 32        |
|                   | for grids in the ParticleBoxArray if dual_grid is true                |             |           |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_grid_size_z   | Maximum number of cells at level 0 in each grid in z-direction        |    Int      | 32        |
|                   | for grids in the ParticleBoxArray if dual_grid is true.               |             |           |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+
