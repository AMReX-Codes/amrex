.. _Chap:CVODE:

CVODE
=====

AMReX supports local ODE integration using the CVODE solver, [1]_ which is part
of the SUNDIALS framework. [2]_ CVODE contains solvers for stiff and non-stiff
ODEs, and as such is well suited for solving e.g., the complex chemistry
networks in combustion simulations, or the nuclear reaction networks in
astrophysical simulations.

Most of CVODE is written in C, but many functions also come with two distinct
Fortran interfaces.  One interface is FCVODE, which is bundled with the stable
release of CVODE.  Its usage is described in the CVODE documentation. [3]_
However, the use of FCVODE is discouraged in AMReX due to its incompatibility
with being used inside OpenMP parallel regions (which is the primary use case
in AMReX applications).

The alternative, and recommended, Fortran interface to  uses the
``iso_c_binding`` feature of the Fortran 2003 standard to implement a direct
interface to the C functions in CVODE.  When compiling CVODE, one need not
build the CVODE library with the FCVODE interface enabled at all.  Rather, the
Fortran 2003 interface to CVODE is provided within AMReX itself.  The
CVODE tutorials provided in AMReX use this new interface.

.. toctree::
   :maxdepth: 1

   CVODE
   SUNDIALS3


.. [1]
   https://computation.llnl.gov/projects/sundials/cvode

.. [2]
   https://computation.llnl.gov/projects/sundials

.. [3]
   https://computation.llnl.gov/sites/default/files/public/cv_guide.pdf
