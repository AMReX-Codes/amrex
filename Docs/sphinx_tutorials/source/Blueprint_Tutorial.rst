.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Blueprint
==========================

These tests, ``AssignMultiLevelDensity`` and ``HeatEquation_EX1_C``,
demonstrate how to convert AMReX Mesh data into an in-memory
Conduit Mesh Blueprint description for consumption by the ALPINE Ascent
in situ visualization and analysis tool.  These are variants, respectively, of 
``amrex/Tests/Particles/AssignMultiLevelDensity`` and
``amrex/Tutorials/Basic/HeatEquation_EX1_C``.

For details about what mesh features are currently supported, see:
amrex/Src/Base/AMReX_Conduit_Blueprint.H

These tests use the interfaces in Src/Base/AMReX_Conduit_Blueprint.H, which
are built when USE_CONDUIT=TRUE. These tests' GNUmakefiles provide a
template of how to enable and link Conduit and Ascent.

For more details about Conduit and Ascent, please see:

Conduit:
  Repo: https://github.com/llnl/conduit
  Docs http://llnl-conduit.readthedocs.io/en/latest/
  Blueprint Docs: http://llnl-conduit.readthedocs.io/en/latest/blueprint.html

Ascent:
  Ascent Repo: http://github.com/alpine-dav/ascent
  Ascent Docs: http://ascent.readthedocs.io/en/latest/

(or ping Cyrus Harrison <cyrush@llnl.gov> or Matt Larsen <<larsen30@llnl.gov>)

