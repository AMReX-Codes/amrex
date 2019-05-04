.. _Chap:Managing the Grid Hierarchy:

Managing the Grid Hierarchy
===========================

There are four separate parts of any strategy to balance computational
work ina  large-scale calculation with hybrid parallelism:

1) Creation of grids -- this includes defining the BoxArray on which
MultiFabs will be built at each level and defining the work estimates

2) Distribution of grids to MPI processes -- this uses the work estimate
defined above (or defaults to work estimate = number of cells) to define
the DistributionMapping with which MultiFabs at that level will be built.

3) When running on multicore machines with OpenMP: creation of grid tiles
(by defining fabarray_mfiter.tile_size), and if relevant, creation of 
particle tiles (by defining particle.tile_size)

4) When running on multicore machines with OpenMP: distribution of tiles to threads

If the application contains task parallelism, load balancing may require special 
care in task scheduling.  See the section on task parallelism.


.. toctree::
   :maxdepth: 1

   ManagingGridHierarchy
