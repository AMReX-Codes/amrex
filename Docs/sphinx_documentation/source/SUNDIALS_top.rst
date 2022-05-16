.. _Chap:SUNDIALS:

SUNDIALS
========

SUNDIALS stands for **SU**\ite of **N**\onlinear and **DI**\fferential/**AL**\gebraic
equation **S**\olvers. It consistes of the following six solvers:

- CVODE, for initial value problems for ODE systems

- CVODES, solves ODE systems and includes sensitivity analysis

- ARKODE, solves initial value ODE problems with Runge-Kutta methods

- IDA, solves initial value problems for differential-algebraic equation systems

- IDAS, solves differential-algebraic equation systems and includes sensitivity analysis

- KINSOL, solves nonlinear algebraic systems


AMReX provides interfaces to the SUNDIALS suite. For time integration, users can
refer to the section :ref:`sec:time_int:sundials` for more information.
In addition, an example code demonstrating time integration with SUNDIALS
can be found in the tutorials at, `SUNDIALS and Time Integrators`_

.. _`SUNDIALS and Time Integrators`: https://amrex-codes.github.io/amrex/tutorials_html/SUNDIALS_Tutorial.html#tutorials-sundials


For more information on SUNDIALS please see
their `readthedocs page <https://sundials.readthedocs.io/en/latest/>`_.
