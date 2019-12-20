.. _developers-repo-structure:

WarpX structure
===============

Repo organization
-----------------

All the WarpX source code is located in ``Source/``. All sub-directories have a pretty straigtforward name. The PIC loop is part of the WarpX class, in function ``WarpX::EvolveEM`` implemented in ``Source/WarpXEvolveEM.cpp``. The core of the PIC loop (i.e., without diagnostics etc.) is in ``WarpX::OneStep_nosub`` (when subcycling is OFF) or ``WarpX::OneStep_sub1`` (when subcycling is ON, with method 1).

Code organization
-----------------

The main WarpX class is WarpX, implemented in ``Source/WarpX.cpp``.

Build system
------------

WarpX uses the AMReX build system. Each sub-folder contains a file ``Make.package`` with the names of header files and source files to include for the build.

WarpX-specific vocabulary
-------------------------

- ``Evolve`` is a generic term to advance a quantity (this comes from AMReX). For instance, ``WarpX::EvolveE(dt)`` advances the electric field for duration ``dt``, ``PhysicalParticleContainer::Evolve(...)`` does field gather + particle push + current deposition for all particles in ``PhysicalParticleContainer``, and ``WarpX::EvolveEM`` is the central ``WarpX`` function that performs 1 PIC iteration.
