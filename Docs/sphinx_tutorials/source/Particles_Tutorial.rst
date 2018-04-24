.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Particles
==========================

There are two tutorials in amrex/Tutorials/Particles that demonstrate the basic usage of 
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

