.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/LinearSolvers
==========================

There are five examples in the Tutorials/LinearSolvers directory. 

``ABecLaplacian_C`` demonstrates how to solve with cell-centered data in a C++ framework.
This example shows how to use either hypre or PETSc as a bottom-solver (or to solve 
the equation at the finest level if you set the "max coarsening level" to 0.

``ABecLaplacian_F`` demonstrates how to solve with cell-centered data using the Fortran interfaces.

``NodalPoisson`` demonstrates how to set up and solve a variable coefficient Poisson equation 
with the rhs and solution data on nodes.

``MAC_Projection_EB`` demonstrates how to set up and perform an EB-aware MAC projection
of a specified velocity field on a staggered grid around nine cylindrical objects.

``Nodal_Projection_EB`` demonstrates how to set up and perform an EB-aware nodal projection
of a specified cell-centered velocity field around nine cylindrical objects.

``MultiComponent`` demonstrates how to solve using a multi-component operator with nodal data.
The operator is of the form :math:`D(\mathbf{\phi})_i = \alpha_{ij}\nabla^2\phi_j`, where :math:`\mathbf{\phi}` is the multi-component solution variable and :math:`\alpha_{ij}` is a matrix of coefficients.
The operator (``MCNodalLinOp``) is implemented using the "reflux-free" method. 
