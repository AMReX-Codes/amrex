.. _Chap:Particles:

Particles
=========

In addition to the tools for working with mesh data described in previous
chapters, AMReX also provides data structures and iterators for performing
data-parallel particle simulations. Our approach is particularly suited to
particles that interact with data defined on a (possibly adaptive)
block-structured hierarchy of meshes. Example applications include
Particle-in-Cell (PIC) simulations, Lagrangian tracers, or particles that exert
drag forces onto a fluid, such as in multi-phase flow calculations. The overall
goals of AMReX's particle tools are to allow users flexibility in specifying
how the particle data is laid out in memory and to handle the parallel
communication of particle data. In the following sections, we
give an overview of AMReX's particle classes and how to use them.


.. toctree::
   :maxdepth: 1

   Particle
