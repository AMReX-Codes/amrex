Parallelization in  WarpX
=========================

When running a simulation, the domain is split into independent 
rectangular sub-domains (called **grids**). This is the way AMReX, a core 
component of WarpX, handles parallelization and/or mesh refinement. Furthermore, 
this decomposition makes load balancing possible: each MPI rank typically computes 
a few grids, and a rank with a lot of work can transfer one or several **grids** 
to their neighbors. 

A user 
does not specify this decomposition explicitly. Instead, the user gives hints to 
the code, and the actual decomposition is determined at runtime, depending on 
the parallelization. The main user-defined parameters are 
``amr.max_grid_size`` and ``amr.blocking_factor``. 

AMReX ``max_grid_size`` and ``blocking_factor``
-----------------------------------------------

* ``amr.max_grid_size`` is the maximum number of points per **grid** along each 
  direction (default ``amr.max_grid_size=32`` in 3D).

* ``amr.blocking_factor``: The size of each **grid** must be divisible by the 
  `blocking_factor` along all dimensions (default ``amr.blocking_factor=8``). 
  Note that the ``max_grid_size`` also has to be divisible by ``blocking_factor``.

These parameters can have a dramatic impact on the code performance. Each 
**grid** in the decomposition is surrounded by guard cells, thus increasing the 
amount of data, computation and communication. Hence having a too small 
``max_grid_size``, may ruin the code performance.

On the other hand, a too-large ``max_grid_size`` is likely to result in a single 
grid per MPI rank, thus preventing load balancing. By setting these two 
parameters, the user wants to give some flexibility to the code while avoiding 
pathological behaviors.

For more information on this decomposition, see the 
`Gridding and Load Balancing <https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html>`__ 
page on AMReX documentation. 

For specific information on the dynamic load balancer used in WarpX, visit the 
`Load Balancing <https://amrex-codes.github.io/amrex/docs_html/LoadBalancing.html>`__ 
page on the AMReX documentation.

The best values for these parameters strongly depends on a number of parameters, 
among which numerical parameters:

* Algorithms used (Maxwell/spectral field solver, filters, order of the 
  particle shape factor)

* Number of guard cells (that depends on the particle shape factor and 
  the type and order of the Maxwell solver, the filters used, `etc.`)

* Number of particles per cell, and the number of species

and MPI decomposition and computer architecture used for the run:

* GPU or CPU

* Number of OpenMP threads

* Amount of high-bandwidth memory.

Below is a list of experience-based parameters 
that were observed to give good performance on given supercomputers.

Rule of thumb for 3D runs on NERSC Cori KNL
-------------------------------------------

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell 
solver on Cori KNL for a well load-balanced problem (in our case laser 
wakefield acceleration simulation in a boosted frame in the quasi-linear 
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=64`` and ``amr.blocking_factor=64`` so that the size of 
  each grid is fixed to ``64**3`` (we are not using load-balancing here).

* **8 MPI ranks per KNL node**, with ``OMP_NUM_THREADS=8`` (that is 64 threads 
  per KNL node, i.e. 1 thread per physical core, and 4 cores left to the 
  system).

* **2 grids per MPI**, *i.e.*, 16 grids per KNL node.
