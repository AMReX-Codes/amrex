.. role:: cpp(code)
   :language: c++

.. _Chap:LinearSolvers:

Linear Solvers
==============

In this Chapter we give an overview of the linear solvers in AMReX.
AMReX supports two solvers where the solution, :math:`\phi`
is either cell-centered or nodal. 
For each case you can solve on a single level (using essentially Dirichlet conditions 
supplied from coarser levels or physical boundary conditions), or you can do a 
multilevel composite solve on a subset of AMR levels (or even all the AMR levels).

The tutorial :cpp:`amrex/Tutorials/LinearSolvers/ABecLaplacian_C` shows how to call the 
cell-centered solver. The tutorial :cpp:`amrex/Tutorials/Basic/HeatEquation_EX3_C` shows 
how to solve the implicit heat equation using the cell-centered solver. Here an 
"ABecLaplacian" has the form 

.. math:: (a_{coeff} A - b_{coeff} \nabla \cdot B \nabla ) \phi = RHS

where :math:`a_{coeff}`, :math:`b_{coeff}` are scalars, :math:`A`, :math:`B`, 
and :math:`RHS` are MultiFabs.
For the cell-centered 
solver, :math:`A` and :math:`RHS` are cell centered, and :math:`B` is an array of 
face-centered MultiFabs since :math:`B \nabla \phi` lives on faces, and the 
divergence lives on cell-centers.

.. toctree::
   :maxdepth: 1

   LinearSolvers
