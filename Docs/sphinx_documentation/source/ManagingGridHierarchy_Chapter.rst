.. role:: cpp(code)
   :language: c++

Gridding and Load Balancing
===========================

AMReX provides a great deal of generality when it comes to how to decompose the 
computational domain into individual logically rectangular grids, and how to distribute
those grids to MPI ranks.  We use the phrase "load balancing" here to refer to the combined process
of grid creation (and re-creation when regridding) and distribution of grids to MPI ranks.

Even for single-level calculations, AMReX provides the flexibility to have different size grids,
more than one grid per MPI rank, and different strategies for distributing the grids to MPI ranks.

For multi-level calculations, the same principles for load balancing apply as in single-level calculations,
but there is additional complexity in how to tag cells for refinement and how to create the 
union of grids at levels > 0 where that union most likely does not cover the computational domain.

See :ref:`sec:grid_creation` for grids are created, i.e. how the :cpp:`BoxArray` on which 
:cpp:`MultiFabs` will be built is defined at each level.

See :ref:`sec:load_balancing` for the strategies AMReX supports for distributing
grids to MPI ranks, i.e. defining the :cpp:`DistributionMapping` with which 
:cpp:`MultiFabs` at that level will be built.  

We also note that we can create separate grids, and map them in different ways to MPI ranks, for 
different types of data in a single calculation.  We refer to this as the "dual grid approach"
and the most common usage is to load balance mesh and particle data separately. See :ref:`sec:dual_grid`
for more about this approach.

When running on multicore machines with OpenMP, we can also control the distribution of 
work by setting the size of grid tiles (by defining :cpp:`fabarray_mfiter.tile_size`), and if relevant, of 
particle tiles (by defining :cpp:`particle.tile_size`).  We can also specify the strategy for assigning 
tiles to OpenMP threads.  See :ref:`sec:basics:mfiter:tiling:` for more about tiling.

.. toctree::
   :maxdepth: 1

   GridCreation
   DualGrid
   LoadBalancing
