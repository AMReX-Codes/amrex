.. role:: cpp(code)
   :language: c++

.. _Chap:LinearSolvers:

Linear Solvers
==============

AMReX supports both single-level solves and composite solves on multiple AMR levels,
with the solution to the linear system defined on cell centers, edges or nodes.
AMReX also supports solution of linear systems with embedded boundaries.
(See chapter :ref:`Chap:EB` for more details on the embedded boundary representation of
complex geometry.) In addition to the iterative solvers discussed in this
chapter, AMReX also supports solving Poisson equations using Fast Fourier
Transform (FFT) (see chapter :ref:`Chap:FFT` for more information).

The default solution technique is geometric multigrid. AMReX also provides
BiCGStab solvers, GMRES, and interfaces to the hypre and PETSc libraries.

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
and :math:`\beta` (which we later refer to as :math:`\sigma`) is cell-centered.

In addition to these solvers, AMReX has support for tensor solves used
to calculate the viscous terms that appear in the compressible Navier-Stokes
equations.  In these solves, all components of the velocity field are solved
for simultaneously.  The tensor solve functionality is only available for
cell-centered velocity. AMReX also supports the curl-curl operator, with
solutions defined on cell edges.

The tutorials in `Linear Solvers`_ show examples of
using the solvers.  The tutorial
``amrex-tutorials/ExampleCodes/Basic/HeatEquation_EX3_C`` shows how to solve the
heat equation implicitly using the solver.  The tutorials also show
how to add linear solvers into the build system.

.. _`Linear Solvers`: https://amrex-codes.github.io/amrex/tutorials_html/LinearSolvers_Tutorial.html

.. toctree::
   :maxdepth: 1

   LinearSolvers
