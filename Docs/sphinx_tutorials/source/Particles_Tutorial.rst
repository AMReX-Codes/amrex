.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Particles
==========================

There are several tutorials in amrex/Tutorials/Particles that demonstrate the basic usage of 
AMReX's particle data structures. 

**ElectrostaticPIC**
--------------------

This tutorial demonstrates how to perform an electrostatic Particle-in-Cell calculation
using AMReX. The code initializes a single particle in a conducting box (i.e. Dirichlet 
zero boundary conditions) that is slightly off-center in one direction. Because of the 
boundary conditions, the particle sees an image charge and is accelerated in this direction.

The code is currently set up to use one level of static mesh refinement. The charge density,
electric field, and electrostatic potential are all defined on the mesh nodes. To solve 
Poisson's equation, we use AMReX's Fortran-based multigrid solver. The Fortran routines for
performing charge deposition, field gathering, and the particle push are all defined in 
:cpp:`electrostatic_pic_2d.f90` and :cpp:`electrostatic_pic_3d.f90` for 2D and 3D, respectively.

The particle container in this example using a Struct-of-Arrays layout, with :cpp:`1 + 2*BL_SPACEDIM`
real components to store the particle weight, velocity, and the electric field interpolated 
to the particle position. To see how to set up such a particle container, see 
:cpp:`ElectrostaticParticleContainer.H`.

**ElectromagneticPIC**
-----------------------

This tutorial shows how to perform an electromagnetic particle-in-cell calculation
using AMReX. Essentially, this is a mini-app version of the WarpX application code.
The electric fields, magnetic fields, and current densities are stored using the
staggered Yee grid, and it solves Maxwell's Equations using the finite-difference
time domain method.

This tutorial also demonstrates how to offload calculations involving particle data
onto the GPU using OpenACC. To compile with GPU support, use the pgi compiler, and set
:cpp:`USE_ACC = TRUE`, and :cpp:`USE_CUDA = TRUE`, :cpp:`USE_OMP = FALSE`. 

You can choose between two problem types by toggling the :cpp:`problem_type` parameter
in the provided inputs file. Choosing the uniform plasma setup provides a nearly
perfectly load balanced problem setup that is useful for performance testing. Choosing
the Langmuir wave problem will automatically compare the simulated fields to the exact
solution.
     
Currently, this tutorial does not use mesh refinement. 
     
**NeighborList**
----------------

This tutorial demonstrates how to have AMReX's particles undergo short-range collisions
with each other. To facilite this, a neighbor list data structure is created, in which
all of the partners that could potentially collide with a given particle are pre-computed.
This is done by first constructing a cell-linked list, and then looping over all 27 neighbor
cells to test for potential collision partners. The Fortran subroutine :cpp:`amrex_compute_forces_nl`
defined in :cpp:`neighbor_list_2d.f90` and :cpp:`neighbor_list_3d.f90` demonstrates how to loop over
the resulting data structure.

The particles in this example store velocity and acceleration in addition to the default
components. They are initially placed at cell centers and given random velocities. When a 
particle reaches the domain boundary, it is specularly reflected back into the domain. To 
see how the particle data structures are set up, see :cpp:`NeighborListParticleContainer.cpp`.

The file called :cpp:`inputs` can be used to run this tutorial with a single level, and
:cpp:`inputs.mr` sets up a run with static mesh refinement.

**CellSortedParticles**
-----------------------

Sometimes, it's useful to sort particles at a finer granularity than grids or tiles. In this
Tutorial, each cell contains a list of particle indices that tell you which particles belong to
that cell. This is useful, for example, in Direct Simulation Monte Carlo calculations, where you want to
potentially interact particles that are in the same cell as each other. Every time the particles move, we check to see
whether it's still in the same cell or not. If it isn't, we mark the particle as unsorted. We then
call Redistribute() as normal, and then insert the unsorted particles into the proper cells. Care is
taken so that, if the Redistribute call changes the order of the particles in the Container, the indices
in the cell lists are updated accordingly. 

This Tutorial is currently single-level only.
