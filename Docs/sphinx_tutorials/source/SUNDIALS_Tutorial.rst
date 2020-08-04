.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/SUNDIALS
==========================

AMReX provides five tutorials in the ``amrex/Tutorials/SUNDIALS`` directory.
There are three versions of ``EX1`` which parallelize differently. ``EX1_C``
packs a box worth of equations into a serial NVector, uses CVODE to solve, and then unpacks
the solution back into the box it came from. ``EX1_CUDA`` uses the cuda NVector implementation
instead. ``EX1_F`` parallelizes over the cells individually. ``EX2_F`` is based on
``fcvRoberts_dns.f`` example code in CVODE. ``EX-CUSOLVER`` uses a Castro-style driver and
tests different ode solving configurations.

See the SUNDIALS_ section of the AMReX documentation for general instructions 
on how to include SUNDIALS in an AMReX application.

.. _SUNDIALS: https://amrex-codes.github.io/amrex/docs_html/SUNDIALS.html
