.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/Amr
=============

For each of these tutorials, 
plotfiles are generated that can be viewed with amrvis2d / amrvis3d
(CCSE's native vis / spreadsheet tool, downloadable separately from ccse.lbl.gov)
or with VisIt.

**Advection_AmrCore**
---------------------

Advection_AmrCore: This tutorial contains an AMR advection code that advects
a single scalar field with a velocity field that is specified on faces.

It is an AMReX based code designed to run in parallel using MPI/OMP.

This example uses source code from the amrex/Src/Base, Boundary, and AmrCore
directories.

Notably, this example does not use source code from amrex/Src/Amr
(see the tutorial Advection_AmrLevel).

The directory Exec/SingleVortex includes a makefile and a sample inputs file.  

**Advection_AmrLevel**
----------------------

Advection_AmrLevel: This tutorial contains an AMR advection code that advects
a single scalar field with a velocity field that is specified on faces.

It is an AMReX based code designed to run in parallel using MPI/OMP.

This example uses source code from the amrex/Src/Base, Boundary, AmrCore, and
Amr directories.

The directories Exec/SingleVortex and Exec/UniformVelocity each include 
a makefile and a sample inputs file.  

**Advection_F**
----------------
This code advects a single scalar field with a velocity
field that is specified on faces.

It is a AMReX based code designed to run in parallel using MPI/OMP.
It uses the Fortran interfaces of AMReX.

The directory Exec/SingleVortex includes a makefile and a sample inputs file.  

**Advection_octree_F**
----------------------

This code advects a single scalar field with a velocity
field that is specified on faces.

It is a AMReX based code designed to run in parallel using MPI/OMP.
It uses the Fortran interfaces of AMReX.  The grids have an octree
structure with a grid size of 8.  No subcycling is used.

The directory Exec/SingleVortex includes a makefile and a sample inputs file.  
