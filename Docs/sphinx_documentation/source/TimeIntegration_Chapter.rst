
.. _sec:basics:timeintegration:

Time Integration
================

AMReX provides a basic explicit time integrator capable of Forward Euler or
both predefined and custom Runge-Kutta schemes designed to advance data on a
particular AMR level by a timestep. This integrator is designed to be flexible,
requiring the user to supply a right-hand side function taking a ``MultiFab``
of state data and filling a ``MultiFab`` of the corresponding right hand side
data. The user simply needs to supply a C++ lambda function to implement
whatever right hand side operations they need.

A Simple Time Integrator Setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is best shown with some sample code that sets up a time integrator and
asks it to step forward by some interval ``dt``. The user needs to supply the
right-hand side function using the ``TimeIntegrator::set_rhs()`` function.

.. highlight:: c++

::

   #include <AMReX_TimeIntegrator.H>

   MultiFab Sold; // MultiFab containing old-time state data
   MultiFab Snew; // MultiFab where we want new-time state data

   // [Fill Sold here]

   // Create a time integrator that will work with
   // MultiFabs with the same BoxArray, DistributionMapping,
   // and number of components as the state_data MultiFab.
   TimeIntegrator<MultiFab> integrator(Sold);

   // Create a function that fills the state BCs and computes the RHS
   auto rhs_fun = [&](MultiFab& rhs, MultiFab& state, const Real time){
       // [Calculate the rhs MultiFab given the state MultiFab]
   };

   // Attach the right hand side function to the integrator
   integrator.set_rhs(source_fun);

   // integrate forward one step from `time` by `dt` to fill Snew
   integrator.advance(Sold, Snew, time, dt);

Picking A Time Integration Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user can customize which integration method they wish to use with a set of
runtime parameters that allow choosing between a simple Forward Euler method or
a generic explicit Runge-Kutta method. If Runge-Kutta is selected, then the user
can choose which of a set of predefined Butcher Tables to use, or can choose to
use a custom table and supply it manually.

The full set of integrator options are detailed as follows:

::

  # INTEGRATION

  ## *** Selecting the integrator backend ***
  ## integration.type can take on the following string or int values:
  ## (without the quotation marks)
  ## "ForwardEuler" or "0" = Native Forward Euler Integrator
  ## "RungeKutta" or "1"   = Native Explicit Runge Kutta
  ## "SUNDIALS" or "2"     = SUNDIALS Integrator
  ## for example:
  integration.type = RungeKutta

  ## *** Parameters Needed For Native Explicit Runge-Kutta ***
  #
  ## integration.rk.type can take the following values:
  ### 0 = User-specified Butcher Tableau
  ### 1 = Forward Euler
  ### 2 = Trapezoid Method
  ### 3 = SSPRK3 Method
  ### 4 = RK4 Method
  integration.rk.type = 3

  ## If using a user-specified Butcher Tableau, then
  ## set nodes, weights, and table entries here:
  #
  ## The Butcher Tableau is read as a flattened,
  ## lower triangular matrix (but including the diagonal)
  ## in row major format.
  integration.rk.weights = 1
  integration.rk.nodes = 0
  integration.rk.tableau = 0.0

.. _sec:time_int:sundials:

Using SUNDIALS
^^^^^^^^^^^^^^

The AMReX Time Integration interface also supports a SUNDIALS backend that
provides explicit, implicit, and implicit-explicit (ImEx) Runge-Kutta methods
as well a multirate (MRI) methods from the ARKODE package in SUNDIALS.
To use SUNDIALS integrators, the user needs to compile AMReX with
``USE_SUNDIALS=TRUE`` and use SUNDIALS v6.0 or later.

The SUNDIALS interface supports ``MultiFab`` or ``Vector<MultiFab>`` data
types. Using a ``Vector<MultiFab>`` permits integrating state data with
different spatial centering (e.g. cell centered, face centered, node centered)
concurrently.

The same code as above can be used with SUNDIALS explicit or implicit
Runge-Kutta methods without any modification. To select a SUNDIALS explicit
Runge-Kutta integrator, one needs only to add the following two input parameters
at runtime:

::

  integration.type = SUNDIALS
  integration.sundials.type = ERK

One can select a different method type by changing the value of
``integration.sundials.type`` to one of the following values:

+------------------------+--------------------------------------------------+
| Input Option           | SUNDIALS Method Type                             |
+========================+==================================================+
| ERK                    | Explicit Runge-Kutta method                      |
+------------------------+--------------------------------------------------+
| DIRK                   | Diagonally Implicit Runge-Kutta method           |
+------------------------+--------------------------------------------------+
| IMEX-RK                | Implicit-Explicit Additive Runge-Kutta method    |
+------------------------+--------------------------------------------------+
| EX-MRI                 | Explicit Multirate Infinitesimal method          |
+------------------------+--------------------------------------------------+
| IM-MRI                 | Implicit Multirate Infinitesimal method          |
+------------------------+--------------------------------------------------+
| IMEX-MRI               | Implicit-Explicit Multirate Infinitesimal method |
+------------------------+--------------------------------------------------+

For ImEx methods, the user needs to supply two right-hand side functions, an
implicit and an explicit function, using the function
``TimeIntegrator::set_imex_rhs()``. Similarly for multirate methods, the user
needs to supply slow and fast right-hand side functions using
``TimeIntegrator::set_rhs()`` to set the slow function and
``TimeIntegrator::set_fast_rhs()`` to set the fast function. With multirate
methods, one also needs to select the fast time scale method type using the
input option ``integration.sundials.fast_type`` which maybe set to ``ERK`` or
``DIRK``.

To select a specific SUNDIALS method use the input option
``integration.sundials.method`` for ERK and DIRK methods as well as the slow
time scale method with an MRI integrator, use ``integration.sundials.method_i``
and ``integration.sundials.method_e`` to set the implicit and explicit method in
an ImEx method, and ``integration.sundials.fast_method`` to set the ERK or DIRK
method used at the fast time scale with an MRI integrator. These options may be
set to any valid SUNDIALS method name, see the following sections in the
SUNDIALS documentation for more details:

* `ERK methods <https://sundials.readthedocs.io/en/latest/arkode/Butcher_link.html#explicit-butcher-tables>`_
* `DIRK methods <https://sundials.readthedocs.io/en/latest/arkode/Butcher_link.html#implicit-butcher-tables>`_
* `ImEx methods <https://sundials.readthedocs.io/en/latest/arkode/Butcher_link.html#additive-butcher-tables>`_
* `MRI methods <https://sundials.readthedocs.io/en/latest/arkode/Usage/MRIStep/MRIStepCoupling.html#mri-coupling-tables>`_

The full set of integrator options are detailed as follows:

::

  # INTEGRATION WITH SUNDIALS

  # *** Select the SUNDIALS integrator backend ***
  integration.type = SUNDIALS

  # *** Select the SUNDIALS method type ***
  # ERK      = Explicit Runge-Kutta method
  # DIRK     = Diagonally Implicit Runge-Kutta method
  # IMEX-RK  = Implicit-Explicit Additive Runge-Kutta method
  # EX-MRI   = Explicit Multirate Infatesimal method
  # IM-MRI   = Implicit Multirate Infatesimal method
  # IMEX-MRI = Implicit-Explicit Multirate Infatesimal method
  integration.sundials.type = ERK

  # *** Select a specific SUNDIALS ERK method ***
  integration.sundials.method = ARKODE_BOGACKI_SHAMPINE_4_2_3

  # *** Select a specific SUNDIALS ImEx method ***
  integration.sundials.method_i = ARKODE_ARK2_DIRK_3_1_2
  integration.sundials.method_e = ARKODE_ARK2_ERK_3_1_2

  # *** Select a specific SUNDIALS MRI method ***
  integration.sundials.method = ARKODE_MIS_KW3
  integration.sundials.fast_method = ARKODE_KNOTH_WOLKE_3_3

The features of this interface evolve with the needs of our codes, so they may
not yet support all SUNDIALS configurations available. If you find you need
SUNDIALS options we have not implemented, please let us know.
