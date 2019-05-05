.. role:: cpp(code)
   :language: c++

Managing the Grid Hierarchy
===========================

Computational load balancing is based on several different steps:

#. The first component of load balancing is the creation of individual grids,
   i.e. defining the :cpp:`BoxArray` on which :cpp:`MultiFabs` will be built at each level

#. The second component is distribution of these grids to MPI processes,
   i.e. defining the :cpp:`DistributionMapping` with which :cpp:`MultiFabs` at that level will be built.
   How we do this depends on what strategy we specify and what work estimate we use

#. When running on multicore machines with OpenMP, we must set the size of grid tiles
   (by defining :cpp:`fabarray_mfiter.tile_size`), and if relevant, of 
   particle tiles (by defining :cpp:`particle.tile_size`).  We must also specify 
   the strategy for assigning tiles to OpenMP threads.

We note that we can load balance the mesh and particle data separately; see :ref:`ss:dual_grid`

.. toctree::
   :maxdepth: 1

   GridCreation
   DualGrid
