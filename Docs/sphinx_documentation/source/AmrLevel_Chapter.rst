.. _Chap:AmrLevel:

Amr Source Code
===============

The source code in ``amrex/Src/Amr`` contains a number of classes, most notably
:cpp:`Amr`, :cpp:`AmrLevel`, and :cpp:`LevelBld`.
These classes provide a more well developed set of tools for writing AMR codes
than the classes created for the Advection_AmrCore tutorial.

-  The :cpp:`Amr` class is derived from :cpp:`AmrCore`, and manages data across the
   entire AMR hierarchy of grids.

-  The :cpp:`AmrLevel` class is a pure virtual class for managing data at a
   single level of refinement.

-  The :cpp:`LevelBld` class is a pure virtual class for defining variable types
   and attributes.

Many of our mature, public application codes contain derived classes that
inherit directly from :cpp:`AmrLevel`. These include:

-  The :cpp:`Castro` class in our compressible astrophysics code, CASTRO,
   (available in the AMReX-Astro/Castro github repository)

-  The :cpp:`Nyx` class in our computational cosmology code, Nyx
   (available in the AMReX-Astro/Nyx github repository).

-  Our incompressible Navier-Stokes code, IAMR
   (available in the AMReX-codes/IAMR github repository) has a pure virtual
   class called :cpp:`NavierStokesBase` that inherits from :cpp:`AmrLevel`, and an
   additional derived class :cpp:`NavierStokes`.

-  Our low Mach number combustion code PeleLM
   (available in the AMReX-Combustion/PeleLM github repository)
   contains a derived class :cpp:`PeleLM` that also inherits from
   :cpp:`NavierStokesBase` (but does not use :cpp:`NavierStokes`).

The tutorial code in ``amrex/Tutorials/Amr/Advection_AmrLevel`` gives a simple
example of a class derived from :cpp:`AmrLevel` that can be used to solve the
advection equation on a subcycling-in-time AMR hierarchy. Note that example is
essentially the same as the ``amrex/Tutorials/Amr/Advection_AmrCore`` tutorial and
documentation in the chapter on :ref:`Chap:AmrCore`, except now we use the
provided libraries in ``amrex/Src/Amr``.

The tutorial code also contains a :cpp:`LevelBldAdv` class (derived from
:cpp:`LevelBld` in the Source/Amr directory). This class is used to define
variable types (how many, nodality, interlevel interpolation stencils, etc.).


.. toctree::
   :maxdepth: 1


   AmrLevel
