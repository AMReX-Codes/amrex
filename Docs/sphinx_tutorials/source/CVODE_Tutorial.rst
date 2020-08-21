.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/CVODE
==========================

There are two CVODE tutorials in the ``amrex/Tutorials/CVODE`` directory, called
``EX1_F`` and ``EX2_F``.  ``EX1_F`` consists of a single ODE that is integrated with
CVODE within each cell of a 3-D grid.  It demonstrates how to initialize the
CVODE solver, how to call the ODE right-hand-side (RHS), and, more importantly,
how to *re-*\ initialize the solver between cells, which avoids allocating and
freeing solver memory between each cell (see the call to ``FCVReInit()`` in the
``integrate_ode.f90`` file in the ``EX1_F`` directory.)

The ``EX2_F`` example demonstrates the slightly more complicated case of
integrating a system of coupled ODEs within each cell.  Similarly to ``EX1_F``,
it provides an RHS and some solver initialization.  However, it also
demonstrates the performance effect of providing an analytic Jacobian matrix
for the system of ODEs, rather than requiring the solver to compute the
Jacobian matrix numerically using a finite-difference approach.  The tutorial
integrates the same system of ODEs on the same 3-D grid, but in one sweep it
instructs CVODE to use the analytic function that computes the Jacobian matrix,
and in the other case, it does not, which requires CVODE to compute it
manually.  One observes a significant performance gain by providing the
analytic Jacobian function.

See the CVODE_ section of the AMReX documentation for general instructions 
on how to include CVODE in an AMReX application.

.. _CVODE: https://amrex-codes.github.io/amrex/docs_html/CVODE.html#id1
