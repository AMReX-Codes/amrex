.. _Chap:AmrLevel:

Amr Source Code
===============

The source code in amrex/Src/Amr contains a number of classes, most notably
Amr, AmrLevel, and LevelBld.
These classes provide a more well developed set of tools for writing AMR codes
than the classes created for the Advection_AmrCore tutorial.

-  The Amr class is derived from AmrCore, and manages data across the
   entire AMR hierarchy of grids.

-  The AmrLevel class is a pure virtual class for managing data at a
   single level of refinement.

-  The LevelBld class is a pure virtual class for defining variable types
   and attributes.

Many of our mature, publicly application codes contain derived classes that inherit directly
from AmrLevel. These include our compressible astrophysics code,
CASTRO,
(available in the AMReX-Astro/Castro github repository)
and our computational cosmology code, Nyx
(available in the AMReX-Astro/Nyx github repository) .
Our incompressible Navier-Stokes code, IAMR
(available in the AMReX-codes/IAMR github repository)
has a pure virtual class called NavierStokesBase that inherits from AmrLevel,
and an additional derived class NavierStokes.
Our low Mach number combustion code PeleLM (not yet public) also inherits
from NavierStokesBase.

The tutorial code in amrex/Tutorials/Amr/Advection_AmrLevel gives a simple
example of a class derived from AmrLevel that can be used to solve
the advection equation on a subcycling-in-time AMR hierarchy. Note that example
is essentially the same as the amrex/Tutorials/Amr/Advection_AmrCore tutorial
and documentation in Chapter `[Chap:AmrCore] <#Chap:AmrCore>`__, except now we use the provided
libraries in Src/Amr.

The tutorial code also contains a LevelBldAdv class (derived from LevelBld in the
Source/Amr directory). This class is used to define variable types (how many, nodality,
interlevel interpolation stencils, etc.).


.. toctree::
   :maxdepth: 1


   AmrLevel
