
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
asks it to step forwards by some interval ``dt``. The user needs to supply at
minimum, the right-hand side function using the ``TimeIntegrator::set_rhs()``
function. By using the ``TimeIntegrator::set_post_update()`` function, a user
can also supply a post update function which is called on state data immediately
before evaluating the right-hand side. This post update function is a good
opportunity to fill boundary conditions for Runge-Kutta stage solution data so that
ghost cells are filled when the right hand side function is called on that solution data.

.. highlight:: c++

::

   #include <AMReX_TimeIntegrator.H>

   MultiFab Sborder; // MultiFab containing old-time state data and ghost cells
   MultiFab Snew;    // MultiFab where we want new-time state data
   Geometry geom;    // The domain (or level) geometry

   // [Fill Sborder here]

   // Create a time integrator that will work with
   // MultiFabs with the same BoxArray, DistributionMapping,
   // and number of components as the state_data MultiFab.
   TimeIntegrator<MultiFab> integrator(Sborder);

   // Create a RHS source function we will integrate
   auto source_fun = [&](MultiFab& rhs, const MultiFab& state, const Real time){
       // User function to calculate the rhs MultiFab given the state MultiFab
       fill_rhs(rhs, state, time);
   };

   // Create a function to call after updating a state
   auto post_update_fun = [&](MultiFab& S_data, const Real time) {
       // Call user function to update state MultiFab, e.g. fill BCs
       post_update(S_data, time, geom);
   };

   // Attach the right hand side and post-update functions
   // to the integrator
   integrator.set_rhs(source_fun);
   integrator.set_post_update(post_update_fun);

   // integrate forward one step from `time` by `dt` to fill S_new
   integrator.advance(Sborder, S_new, time, dt);


Picking A Time Integration Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user can customize which integration method they wish to use with a set of
runtime parameters that allow choosing between a simple Forward Euler method or
a generic explicit Runge-Kutta method. If Runge-Kutta is selected, then the user
can choose which of a set of predefined Butcher Tables to use, or can choose to
use a custom table and supply it manually.

When AMReX is compiled with SUNDIALS v.6 or later, the user also has an option
to use the SUNDIALS ARKODE integrator as a backend for the AMReX Time Integrator
class. The features of this interface evolve with the needs of our codes, so
they may not yet support all SUNDIALS configurations available. If you find you
need SUNDIALS options we have not implemented, please let us know.

The options are detailed as follows:

::

  # INTEGRATION

  ## *** Selecting the integrator backend ***
  ## integration.type can take on the following string or int values:
  ## (without the quotation marks)
  ## "ForwardEuler" or "0" = Native Forward Euler Integrator
  ## "RungeKutta" or "1"   = Native Explicit Runge Kutta
  ## "SUNDIALS" or "2"     = SUNDIALS ARKODE Integrator
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

  ## *** Parameters Needed For SUNDIALS ARKODE Integrator ***
  ## integration.sundials.strategy specifies which ARKODE strategy to use.
  ## The available options are (without the quoatations):
  ## "ERK" = Explicit Runge Kutta
  ## "MRI" = Multirate Integrator
  ## "MRITEST" = Tests the Multirate Integrator by setting a zero-valued fast RHS function
  ## for example:
  integration.sundials.strategy = ERK
