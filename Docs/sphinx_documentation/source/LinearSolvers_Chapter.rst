.. role:: cpp(code)
   :language: c++

.. _Chap:LinearSolvers:

Linear Solvers
==============

AMReX supports both single-level solves and composite solves on multiple AMR levels,
with the scalar solution to the linear system defined on either cell centers or nodes.
AMReX also supports solution of linear systems with embedded boundaries. 
(See chapter :ref:`Chap:EB` for more details on the embedded boundary representation of
complex geometry.)

The default solution technique is geometric multigrid, but AMReX includes native
BiCGStab solvers for a single level as well as interfaces to the hypre library.

In this Chapter we give an overview of the linear solvers in AMReX
that solve linear systems in the canonical form

.. math:: (A \alpha - B \nabla \cdot \beta \nabla ) \phi = f,
  :label: eqn::abeclap

where :math:`A` and :math:`B` are scalar constants, 
:math:`\alpha` and :math:`\beta` are scalar fields,
:math:`\phi` is the unknown, 
and :math:`f` is the right-hand side of the equation.  Note
that Poisson's equation :math:`\nabla^2 \phi = f` is a special case
of the canonical form.  The solution :math:`\phi` is at either
cell centers or nodes.  

For the cell-centered solver, :math:`\alpha`, :math:`\phi` and :math:`f` 
are represented by cell-centered MultiFabs,
and :math:`\beta` is represented by ``AMREX_SPACEDIM`` face type
MultiFabs, i.e.  there are separate MultiFabs for the :math:`\beta` 
coefficient in each coordinate direction.

For the nodal solver, :math:`A` and :math:`\alpha` are assumed to be zero, 
:math:`\phi` and :math:`f` are nodal,
and :math:`\beta` (whcih we later refer to as :math:`\sigma`) is cell-centered.  

In addition to these solvers, AMReX has support for tensor solves used
to calculate the viscous terms that appear in the compressible Navier-Stokes
equations.  In these solves, all components of the velocity field are solved
for simultaneously.  The tensor solve functionality is only available for
cell-centered velocity.

The tutorials in ``amrex/Tutorials/LinearSolvers/`` show examples of
using the solvers.  The tutorial
``amrex/Tutorials/Basic/HeatEquation_EX3_C`` shows how to solve the
heat equation implicitly using the solver.  The tutorials also show
how to add linear solvers into the build system.

.. toctree::
   :maxdepth: 1

   LinearSolvers
