.. role:: cpp(code)

.. _Chap:AmrLevel:

Amr Source Code
===============

The source code in ``amrex/Src/Amr`` contains a number of classes, most notably
:cpp:`Amr`, :cpp:`AmrLevel`, and :cpp:`LevelBld`.
These classes provide a more well-developed set of tools for writing AMR codes
than the classes created for the AdvectionAmrCore tutorial.

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

The tutorial code in ``amrex-tutorials/AMR/WaveAmrLevel`` gives a simple
example of a class derived from :cpp:`AmrLevel` that solves the wave
equation on a subcycling-in-time AMR hierarchy, and
``amrex-tutorials/Applications/CNS`` is a complete compressible flow code
built on these classes. Compare them with the `AdvectionAmrCore`_ tutorial
and the chapter on :ref:`Chap:AmrCore`, which build the same kind of
simulation directly on :cpp:`AmrCore` without the classes in ``amrex/Src/Amr``.

.. _`AdvectionAmrCore`: https://amrex-codes.github.io/amrex/tutorials_html/AMR.html#amr-advectionamrcore

The tutorial codes also contain a class derived from :cpp:`LevelBld` (in
``amrex/Src/Amr``). This class is used to define
variable types (how many, nodality, interlevel interpolation stencils, etc.).


.. toctree::
   :maxdepth: 1


   AmrLevel
