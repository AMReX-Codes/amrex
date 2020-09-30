.. _Chap:SUNDIALS:

SUNDIALS
========

AMReX supports local ODE integration using the ARKode [1]_ and CVODE [2]_
time integrators which are part of the SUNDIALS framework [3]_. ARKode
and CVODE contains solvers for stiff and non-stiff ODEs, and as such they
are well suited for solving e.g., the complex chemistry networks in combustion
simulations, or the nuclear reaction networks in astrophysical simulations.

Most of SUNDIALS is written in C, but it is distributed with Fortran
interfaces that use the ``iso_c_binding`` feature of the Fortran 2003 standard.
AMReX supports these Fortran 2003 interfaces and they are used in the AMReX
SUNDIALS 5 tutorials.

AMReX currently supports SUNDIALS version 5, and for CVODE only, a legacy
interface to SUNDIALS 2.7 which is the version available in the ``cray-tpsl``
system module made available on Cray systems.


.. toctree::
   :maxdepth: 2

   SUNDIALS
   SUNDIALS_CVODE

.. [1]
   https://computation.llnl.gov/projects/sundials/arkode

.. [2]
   https://computation.llnl.gov/projects/sundials/cvode

.. [3]
   https://computation.llnl.gov/projects/sundials

