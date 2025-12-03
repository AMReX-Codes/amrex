.. _Chap:Particles:

Particles
=========

In addition to the tools for working with mesh data described in previous chapters,
AMReX provides data structures and iterators for performing data-parallel particle simulations.
While these tool can be used to implement pure particle methods, they are particularly
suited for methods in which particles interact with data defined on mesh or
hierarchy of meshes. Example applications include Particle-in-Cell (PIC) simulations,
Lagrangian tracers, and particles that exert drag forces onto a fluid. The overall
goals of AMReX's particle tools are to allow users to express a variety of useful
operations on particle data, including particle-mesh and particle-particle operations,
in a performance-portable and scalable manner. In the following sections, we
give an overview of AMReX's particle classes and how to use them.


.. toctree::
   :maxdepth: 1

   Particle
