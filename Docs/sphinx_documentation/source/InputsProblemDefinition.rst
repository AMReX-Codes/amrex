.. _sec:inputs:pd:

Problem Definition
==================

The following inputs must be preceded by "amr."

+-------------------+-----------------------------------------------------------------------+-------------+-----------+
|                   | Description                                                           |   Type      | Default   |
+===================+=======================================================================+=============+===========+
| n_cell            | Number of cells at level 0 in each coordinate direction               | Int Int Int | None      |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+
| max_level         | Maximum level of refinement allowed (0 when single-level)             |    Int      | None      |
+-------------------+-----------------------------------------------------------------------+-------------+-----------+

The following inputs must be preceded by "geometry."

+-----------------+-----------------------------------------------------------------------+-------------+-----------+
|                 | Description                                                           |   Type      | Default   |
+=================+=======================================================================+=============+===========+
| coord_sys       | 0 for Cartesian                                                       |   Int       |   0       |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| is_periodic     | 1 for true, 0 for false (one value for each coordinate direction)     |   Ints      | 0 0 0     |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_lo         | Low corner of physical domain (physical not index space)              |   Reals     | 0 0 0     |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_hi         | High corner of physical domain (physical not index space)             |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+
| prob_extent     | Extent of physical domain (physical not index space)                  |   Reals     | None      |
+-----------------+-----------------------------------------------------------------------+-------------+-----------+

Note that internally ``prob_lo`` and ``prob_hi`` are the variables carried by the ``Geometry`` class.
In the inputs file (or command line), one can specify
1) ``geometry.prob_hi`` only or
2) ``geometry.prob_extent`` only or
3) ``geometry.prob_lo`` and ``geometry.prob_hi`` or
4) ``geometry.prob_lo`` and ``geometry.prob_extent``.
If ``geometry.prob_lo`` is not specified then it will be 0 in each coordinate direction.
If ``geometry.prob_extent`` is specified (and ``geometry.prob_hi`` is not) then internally
"prob_hi" will be set to "prob_lo" + "prob_extent".

