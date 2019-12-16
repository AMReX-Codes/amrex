.. WarpX documentation master file, created by
   sphinx-quickstart on Tue Mar  7 22:08:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

WarpX documentation
===================

.. warning::

   This is an **alpha release** of WarpX.
   The code is still in active development.
   Robustness and performance may fluctuate at this stage.
   The input and output formats may evolve.

WarpX is an advanced **electromagnetic Particle-In-Cell** code.

It supports many features including:

    - Perfectly-Matched Layers (PML)
    - Boosted-frame simulations
    - Mesh refinement

For details on the algorithms that WarpX implements, see the section :doc:`theory/theory`.

In addition, WarpX is a highly-parallel and highly-optimized code
and features hybrid OpenMP/MPI parallelization, advanced vectorization
techniques and load balancing capabilities.

In order to learn to use the code, please see the sections below:

.. toctree::
   :maxdepth: 1

   building/building
   running_cpp/running_cpp
   running_python/running_python
   visualization/visualization
   theory/theory
   developers/developers
   api/library_root
