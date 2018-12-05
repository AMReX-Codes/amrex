.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/LinearSolvers
==========================

There are three examples in the Tutorials/LinearSolvers directory. 

``ABecLaplacian_C`` demonstrates how to solve with cell-centered data in a C++ framework.
This example shows how to use either hypre or PETSc as a bottom-solver (or to solve 
the equation at the finest level if you set the "max coarsening level" to 0.

``ABecLaplacian_F`` demonstrates how to solve with cell-centered data using the Fortran interfaces.

``NodalPoisson`` demonstrates how to solve with nodal data using the C++ framework.

