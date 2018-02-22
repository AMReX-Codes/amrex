.. _Chap:AmrCore:

AmrCore Source Code
===================

In this Chapter we give an overview of functionality contained in the
``amrex/Src/AmrCore`` source code.  This directory contains source code for the
following:

-  Storing information about the grid layout and processor distribution mapping
   at each level of refinement.

-  Functions to create grids at different levels of refinement, including
   tagging operations.

-  Operations on data at different levels of refinement, such as interpolation
   and restriction operators.

-  Flux registers used to store and manipulate fluxes at coarse-fine
   interfaces.

-  Particle support for AMR (see :ref:`Chap:Particles`).

There is another source directory, ``amrex/Src/Amr/``, which contains
additional classes used to manage the time-stepping for AMR simulations.
However, it is possible to build a fully adaptive, subcycling-in-time
simulation code without these additional classes.

In this Chapter, we restrict our use to the ``amrex/Src/AmrCore`` source code
and present a tutorial that performs an adaptive, subcycling-in-time simulation
of the advection equation for a passively advected scalar.  The accompanying
tutorial code is available in ``amrex/Tutorials/Amr/Advection_AmrCore`` with
build/run directory ``Exec/SingleVortex``. In this example, the velocity field
is a specified function of space and time, such that an initial Gaussian
profile is displaced but returns to its original configuration at the final
time.  The boundary conditions are periodic and we use a refinement ratio of
:math:`r=2` between each AMR level. The results of the simulation in
two-dimensions are depicted in the Table showing the :ref:`SingleVortex
Tutorial<fig:Adv>`.




.. toctree::
   :maxdepth: 1

   AmrCore
