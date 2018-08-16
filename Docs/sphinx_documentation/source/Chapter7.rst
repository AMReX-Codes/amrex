.. role:: cpp(code)
   :language: c++

.. _Chap:LinearSolvers:

Linear Solvers
==============

In this Chapter we give an overview of the linear solvers in AMReX
that solve linear systems in the canonical form

.. math:: (\alpha A - \beta \nabla \cdot B \nabla ) \phi = f,
  :label: eqn::abeclap

where :math:`\alpha` and :math:`\beta` are scalar constants,
:math:`A` and :math:`B` are scalar fields, :math:`\phi` is the
unknown, and :math:`f` is the right-hand side of the equation.  Note
that Poisson's equation :math:`\nabla^2 \phi = f` is a special case
the canonical form.  The sought solution :math:`\phi` is at either
cell centers or nodes.  For the cell-centered solver, :math:`A`,
:math:`\phi` and :math:`f` are represented by cell-centered MultiFabs,
and :math:`B` is represented by ``AMREX_SPACEDIM`` face type
MultiFabs.  That is there are 1, 2 and 3 MultiFabs representing
:math:`B` for 1, 2 and 3D, respectively.  For the nodal solver,
:math:`A` is assumed to be zero, :math:`\phi` and :math:`f` are nodal,
and :math:`B` is cell-centered.  AMReX supports both single-level
solve and composite solve on multiple AMR levels.

The tutorials in ``amrex/Tutorials/LinearSolvers/`` show examples of
using the solvers.  The tutorial
``amrex/Tutorials/Basic/HeatEquation_EX3_C`` shows how to solve the
heat equation implicitly using the solver.  The tutorials also show
how to add linear solvers into the build system.

.. toctree::
   :maxdepth: 1

   LinearSolvers
