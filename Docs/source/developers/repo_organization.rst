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

WarpX uses the AMReX build system (GNUMake).
Each sub-folder contains a file ``Make.package`` with the names of source files (``.cpp``) that are added to the build.
Do not list header files (``.H``) here.

C++ Includes
------------

All WarpX header files need to be specified relative to the ``Source/`` directory.

- e.g. ``#include "Utils/WarpXConst.H"``
- files in the same directory as the including header-file can be included with ``#include "FileName.H"``

The `include order <https://github.com/ECP-WarpX/WarpX/pull/874#issuecomment-607038803>`_ and `proper quotation marks <https://gcc.gnu.org/onlinedocs/cpp/Include-Syntax.html>`_ are:

1. In a ``<MyName>.cpp`` file: ``#include "<MyName>.H"`` (its header) then
2. (further) WarpX header files ``#include "..."`` then
3. PICSAR and AMReX header files ``#include <...>`` then
4. other third party includes ``#include <...>`` then
5. standard library includes, e.g. ``#include <vector>``

For details why this is needed, please see `PR #874 <https://github.com/ECP-WarpX/WarpX/pull/874#issuecomment-607038803>`_, the `LLVM guidelines <https://llvm.org/docs/CodingStandards.html#include-style>`_, and `include-what-you-use <https://github.com/include-what-you-use/include-what-you-use/blob/master/docs/WhyIWYU.md>`_.

WarpX-specific vocabulary
-------------------------

- ``Evolve`` is a generic term to advance a quantity (this comes from AMReX). For instance, ``WarpX::EvolveE(dt)`` advances the electric field for duration ``dt``, ``PhysicalParticleContainer::Evolve(...)`` does field gather + particle push + current deposition for all particles in ``PhysicalParticleContainer``, and ``WarpX::EvolveEM`` is the central ``WarpX`` function that performs 1 PIC iteration.
