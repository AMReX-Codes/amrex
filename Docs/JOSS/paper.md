---
title: 'AMReX: a framework for block-structured adaptive mesh refinement (AMR)'

tags:
- C++
- adaptive mesh refinement 
- block-structured
- finite volume
- partial differential equations

authors:
- name: Weiqun Zhang
  orcid: 0000-0001-8092-1974
  affiliation: 1
- name: Ann Almgren
  orcid: 0000-0003-2103-312X
  affiliation: 1
- name: Vince Beckner
  affiliation: 1
- name: Johannes Blaschke
  orcid: 0000-0002-6024-3990
  affiliation: 1
- name: Cy Chan
  orcid: 0000-0001-6881-827X
  affiliation: 2
- name: Marcus Day
  orcid: 0000-0002-1711-3963
  affiliation: 1
- name: Brian Friesen
  orcid: 0000-0002-1752-1631
  affiliation: 3
- name: Kevin Gott
  orcid: 0000-0003-3244-5525
  affiliation: 3
- name: Daniel Graves
  orcid: 0000-0001-9730-7217
  affiliation: 2
- name: Maximilian Katz
  orcid: 0000-0003-0439-4556
  affiliation: 4
- name: Andrew Myers
  orcid: 0000-0001-8427-8330
  affiliation: 1
- name: Tan Nguyen
  orcid: 0000-0003-3748-403X
  affiliation: 2
- name: Andrew Nonaka
  orcid: 0000-0003-1791-0265
  affiliation: 1
- name: Michele Rosso
  orcid: 0000-0001-8126-7425
  affiliation: 1
- name: Samuel Williams
  orcid: 0000-0002-8327-5717
  affiliation: 2
- name: Michael Zingale
  orcid: 0000-0001-8401-030X
  affiliation: 5

affiliations:
- name: Center for Computational Sciences and Engineering (CCSE), Lawrence Berkeley National Laboratory
  index: 1
- name: Computational Research Division, Lawrence Berkeley National Laboratory
  index: 2
- name: National Energy Research Scientific Computing Center (NERSC)
  index: 3
- name: NVIDIA Corporation
  index: 4
- name: Department of Physics and Astronomy, Stony Brook University
  index: 5

date: 22 March 2019

bibliography: paper.bib
---

# Summary

`AMReX` is a C++ software framework that supports the development of 
block-structured adaptive mesh refinement (AMR) algorithms for solving 
systems of partial differential equations (PDEs) with complex boundary 
conditions on current and emerging architectures. 
AMR reduces the computational cost and memory footprint
compared to a uniform mesh while preserving the essentially 
local descriptions of different physical processes in complex multiphysics algorithms. 
The origins of AMReX trace back to the BoxLib software framework.

AMReX supplies data containers and iterators for mesh-based fields,
particle data and irregular embedded boundary (cut cell) representations of complex geometries.
The mesh-based data can be defined on cell centers, cell faces or cell corners (nodes).
Coordinate systems include 1D Cartesian or spherical; 2D Cartesian or cylindrical (r-z); and 3D Cartesian.
Both particles and embedded boundary representations introduce additional irregularity and complexity to the
way data is stored and operated on, requiring special attention in the presence of the dynamically
changing hierarchical mesh structure and AMR timestepping approaches.

AMReX supports a number of different time-stepping strategies
and spatial discretizations.  Solution strategies supported by AMReX range
from level-by-level approaches (with or without subcycling in time) with
multilevel synchronization to full-hierarchy approaches, and any combination thereof.

Block-structured AMR provides the basis for the temporal and spatial 
discretization strategy for a large number of applications.
Current AMReX-based application codes target accelerator design, astrophysics, 
combustion, cosmology, microfluidics, materials science and multiphase flow. 

AMReX is part of the 2018 xSDK software release.

# Particles

AMReX provides data structures and iterators for performing data-parallel particle simulations. 
Our approach is particularly suited to particles that interact with data defined on a (possibly adaptive) 
block-structured hierarchy of meshes. Example applications include Particle-in-Cell (PIC) simulations, 
Lagrangian tracers, or particles that exert drag forces onto a fluid, such as in multi-phase flow calculations. 
The overall goals of AMReX’s particle tools are to allow users flexibility in specifying how the particle data 
is laid out in memory and to handle the parallel communication of particle data. 

# Linear Solvers

AMReX includes native linear solvers for parabolic and elliptic equations.  Solution procedures
include geometric multigrid and BiCGStab iterative solvers; interfaces to external hypre and
PETSc solvers are also provided.

# I/O and Post-processing

AMReX has native I/O for checkpointing and for reading and writing plotfiles for post-processing
analysis or visualization.   AMReX also supplies interfaces to HDF5.  The AMReX plotfile format
is supported by VisIt, Paraview, and yt.   AMReX also has linkages to external routines through
both Conduit and SENSEI.

# Asynchronous Iterators

Hiding communication overheads via overlapping communication with computation requires 
a sufficiently large amount of task parallelism. This problem is even more challenging due to 
various types of tasks in an AMReX program, including data parallel tasks 
(same workload on different data partitions) and control parallel tasks 
(different types of workload). AMReX’s asynchronous iterators can facilitate the 
job of identifying tasks in the applications using two types of iterators.

# Fork-Join Support

An AMReX-based application consists of a set of MPI ranks cooperating together on distributed data. 
Typically, all of the ranks in a job compute in a bulk-synchronous, data-parallel fashion, 
where every rank does the same sequence of operations, each on different parts of the distributed data.

The fork-join functionality in AMReX allows the user to divide the job’s MPI ranks into subgroups (i.e. fork) 
and assign each subgroup an independent task to compute in parallel with each other. 
After all of the forked child tasks complete, they synchronize (i.e. join), and the parent task continues execution as before.

The fork-join operation can also be invoked in a nested fashion, creating a hierarchy of fork-join operations, 
where each fork further subdivides the ranks of a task into child tasks. 
This approach enables heterogeneous computation and reduces the strong scaling penalty for 
operations with less inherent parallelism or with large communication overheads.

# Fortran interfaces

The core of AMReX is written in C++. For Fortran users who want to write all of their 
programs in Fortran, AMReX provides Fortran interfaces around much of AMReX's core functionality.

# Parallelism

AMReX’s GPU strategy focuses on providing performant GPU support with 
minimal changes to AMReX-based application codes and maximum flexibility. 
This allows application teams to get running on GPUs quickly while allowing 
long term perfomance tuning and programming model selection. 
AMReX currently uses CUDA for GPUs, but application teams can use CUDA, CUDA Fortran, 
OpenACC or OpenMP in their individual codes.  AMReX will support non-CUDA strategies 
as appropriate.

When running on CPUs, AMReX uses an MPI+X strategy where the X threads are used to perform 
parallelization techniques like tiling. The most common X is OpenMP. On GPUs, AMReX requires CUDA 
and can be further combined with other parallel GPU languages, including OpenACC and OpenMP, 
to control the offloading of subroutines to the GPU. This MPI+CUDA+X GPU strategy has been developed 
to give users the maximum flexibility to find the best combination of portability, 
readability and performance for their applications

# Tutorials

Demonstrations of using AMReX for building application codes are provided in the AMReX Tutorials section.
Examples include a Particle-in-Cell (PIC) code, a compressible Navier-Stokes solver in complex geometry, 
advection-diffusion solvers,  support for spectral deferred corrections time-stepping, and much more.

# AMReX Profiling Tools

AMReX-based application codes can be instrumented using AMReX-specific performance profiling tools that take 
into account the hierarchical nature of the mesh in most AMReX-based applications. 
These codes can be instrumented for varying levels of profiling detail.

# Acknowledgements

The development of AMReX was supported by the
U.S. Department of Energy, Office of Science, 
Office of Advanced Scientific Computing Research, 
Applied Mathematics program under contract number DE-AC02005CH11231,
and by the  
Exascale Computing Project (17-SC-20-SC), a collaborative effort of the 
U.S. Department of Energy Office of Science and the National Nuclear Security Administration.

# References
