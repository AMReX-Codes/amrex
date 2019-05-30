.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:dual_grid:

Dual Grid Approach
------------------

In AMReX-based applications that have both mesh data and particle data,
the mesh work and particle work have very different requirements for load balancing.

Rather than using a combined work estimate to create the same grids for mesh and particle
data, we have the option to pursue a "dual grid" approach.

With this approach the mesh (:cpp:`MultiFab`) and particle (:cpp:`ParticleContainer`) data 
are allocated on different :cpp:`BoxArrays` with different :cpp:`DistributionMappings`.

This enables separate load balancing strategies to be used for the mesh and particle work.

The cost of this strategy, of course, is the need to copy mesh data onto temporary 
:cpp:`MultiFabs` defined on the particle :cpp:`BoxArrays` when mesh-particle communication
is required.

