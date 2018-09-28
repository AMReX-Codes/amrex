.. WarpX documentation master file, created by
   sphinx-quickstart on Tue Mar  7 22:08:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

WarpX documentation
===================

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

.. warning::

   WarpX is currently in active development, and this release
   is an **alpha release**.
   The input and output formats, as well as performances, 
   are not guaranteed to be stable at this stage.



.. toctree::
   :maxdepth: 1

   installation
   running_cpp/running_cpp
   running_python/running_python
   visualization/visualization
   theory/theory
