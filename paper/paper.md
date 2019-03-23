---
title: 'AMReX: a framework for block-structured adaptive mesh refinement'

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
- name: John Bell
  orcid: 0000-0002-5749-334X
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
  orcid: 0000-0002-1572-1631
  affiliation: 3
- name: Kevin Gott
  orcid: 0000-0003-3244-5525
  affiliation: 3
- name: Daniel Graves
  orcid: 0000-0001-9730-7217
  affiliation: 2
- name: Max P. Katz
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

AMReX is a C++ software framework that supports the development of 
block-structured adaptive mesh refinement (AMR) algorithms for solving 
systems of partial differential equations (PDEs) with complex boundary 
conditions on current and emerging architectures.  

Block-structured AMR provides the basis for the temporal and spatial 
discretization strategy for a large number of applications.
AMR reduces the computational cost 
and memory footprint compared to a uniform mesh while preserving the
local descriptions of different physical processes in complex multiphysics algorithms. 
Current AMReX-based application codes span a number of areas; in particular 
the AMReX-Astro GitHub repository holds a number of astrophysical modeling tools based on AMReX [@Zingale_2018].
The origins of AMReX trace back to the BoxLib [@BoxLib] software framework.

AMReX supports a number of different time-stepping strategies
and spatial discretizations.  Solution strategies supported by AMReX range
from level-by-level approaches (with or without subcycling in time) with
multilevel synchronization to full-hierarchy approaches, and any combination thereof.
User-defined kernels that operate 
on patches of data can be written in C++ or Fortran; there is also a Fortran-interface functionality
which wraps the core C++ data structures and operations in Fortran wrappers so that an application
code based on AMReX can be written entirely in Fortran.

AMReX developers believe that interoperability is an important feature of sustainable software.
AMReX has examples of interfaces to other popular software packages such as SUNDIALS,
PETSc and hypre, and is part of the 2018 xSDK software release thus installable with Spack.  

### Mesh and Particle Data

AMReX supplies data containers and iterators for mesh-based fields and particle data.
The mesh-based data can be defined on cell centers, cell faces or cell corners (nodes).
Coordinate systems include 1D Cartesian or spherical; 2D Cartesian or cylindrical (r-z); and 3D Cartesian.

AMReX provides data structures and iterators for performing data-parallel particle simulations. 
The approach is particularly suited to particles that interact with data defined on a (possibly adaptive) 
block-structured hierarchy of meshes. Example applications include those that use Particle-in-Cell (PIC) methods, 
Lagrangian tracers, or solid particles that exchange momentum with the surrounding fluid through drag forces.
AMReX’s particle implementation allows users flexibility in specifying how the particle data 
is laid out in memory and in choosing how to optimize parallel communication of particle data. 

### Complex Geometries

AMReX provides support for discretizing complex geometries using the cut cell / embedded boundary 
approach.   This requires additional data structures for holding face apertures and normals as well
as volume fractions.   Support for operations on the mesh hierarchy including cut cells is enabled 
through the use of specialized discretizations at and near cut cells, and masks to ensure that only
values in the valid domain are computed.  Examples are provided in the tutorials.

### Parallelism

AMReX’s GPU strategy focuses on providing performant GPU support with 
minimal changes to AMReX-based application codes and maximum flexibility. 
This allows application teams to get running on GPUs quickly while allowing 
long term performance tuning and programming model selection. 
AMReX currently uses CUDA for GPUs, but application teams can use CUDA, CUDA Fortran, 
OpenACC or OpenMP in their individual codes.  AMReX will support non-CUDA strategies 
as appropriate.

When running on CPUs, AMReX uses an MPI+X strategy where the X threads are used to perform 
parallelization techniques like tiling. The most common X is OpenMP. On GPUs, AMReX requires CUDA 
and can be further combined with other parallel GPU languages, including OpenACC and OpenMP, 
to control the offloading of subroutines to the GPU. This MPI+CUDA+X GPU strategy has been developed 
to give users the maximum flexibility to find the best combination of portability, 
readability and performance for their applications.

### Asynchronous Iterators and Fork-Join Support

AMReX includes a runtime system that can execute asynchronous AMReX-based applications efficiently
on large-scale systems.  The runtime system constructs a task dependency graph for the 
whole coarse time step and executes it asynchronously to the completion of the step.  There
is also support for more user-specific algorithms such as asynchronous filling of ghost cells
across multiple ranks, including interpolation of data in space and time.

In addition, AMReX has support for fork-join functionality. During a run of an AMReX-based application,
the user can divide the MPI ranks into subgroups (i.e. fork) and assign each subgroup an independent task
to compute in parallel with each other.
After all of the forked child tasks complete, they synchronize (i.e. join), and the parent task continues execution as before.
The fork-join operation can also be invoked in a nested fashion, creating a hierarchy of fork-join operations, 
where each fork further subdivides the ranks of a task into child tasks. 
This approach enables heterogeneous computation and reduces the strong scaling penalty for 
operations with less inherent parallelism or with large communication overheads.

### Linear Solvers

AMReX includes native linear solvers for parabolic and elliptic equations.  Solution procedures
include geometric multigrid and BiCGStab iterative solvers; interfaces to external hypre and
PETSc solvers are also provided.   The linear solvers operate on regular mesh data as well 
as data with cut cells.

### I/O and Post-processing

AMReX has native I/O for checkpointing and for reading and writing plotfiles for post-processing
analysis or visualization.   AMReX also supplies interfaces to HDF5.  The AMReX plotfile format
is supported by VisIt, Paraview, and yt.   AMReX also has linkages to external routines through
both Conduit and SENSEI.

### Documentation, Tutorials and Profiling Tools

Extensive documentation of core AMReX functionality is available online, and many of the application
codes based on AMReX are publicly available as well.  Smaller examples of using AMReX for building application codes 
are provided in the AMReX Tutorials section.
Examples include a Particle-in-Cell (PIC) code, a compressible Navier-Stokes solver in complex geometry, 
advection-diffusion solvers,  support for spectral deferred corrections time-stepping, and much more.

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
