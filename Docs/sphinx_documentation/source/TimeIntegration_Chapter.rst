
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

.. _sec:time_int:sundials:

Using SUNDIALS
^^^^^^^^^^^^^^

The AMReX Time Integration interface also supports a SUNDIALS backend that
wraps both the explicit Runge-Kutta (ERK) and multirate (MRI) integration
schemes in the SUNDIALS ARKODE package. To use either of them, the user needs
to compile AMReX with `USE_SUNDIALS=TRUE` and use SUNDIALS v. 6.0 or later.

There are only minor changes to the code above required to use the SUNDIALS
interface. The first change is that the integration datatype is now a
`Vector<MultiFab>` type instead of simply `MultiFab`. The reason for
introducing a `Vector<MultiFab>` in this case, is to permit integrating state
data with different spatial centering (e.g. cell centered, face centered, node
centered) concurrently. Shown here is sample code equivalent to the code above,
suitable for the SUNDIALS explicit Runge-Kutta integrator:

.. highlight:: c++

::

   #include <AMReX_TimeIntegrator.H>

   Vector<MultiFab> Sborder; // MultiFab(s) containing old-time state data and ghost cells
   Vector<MultiFab> Snew;    // MultiFab(s) where we want new-time state data
   Geometry geom;    // The domain (or level) geometry

   // [Fill Sborder here]

   // Create a time integrator that will work with
   // MultiFabs with the same BoxArray, DistributionMapping,
   // and number of components as the state_data MultiFab.
   TimeIntegrator<Vector<MultiFab> > integrator(Sborder);

   // Create a RHS source function we will integrate
   auto source_fun = [&](Vector<MultiFab>& rhs, const Vector<MultiFab>& state, const Real time){
       // User function to calculate the rhs MultiFab given the state MultiFab
       fill_rhs(rhs, state, time);
   };

   // Create a function to call after updating a state
   auto post_update_fun = [&](Vector<MultiFab>& S_data, const Real time) {
       // Call user function to update state MultiFab, e.g. fill BCs
       post_update(S_data, time, geom);
   };

   // Attach the right hand side and post-update functions
   // to the integrator
   integrator.set_rhs(source_fun);
   integrator.set_post_update(post_update_fun);

   // integrate forward one step from `time` by `dt` to fill S_new
   integrator.advance(Sborder, S_new, time, dt);

Afterwards, to select the ERK integrator, one needs only to add the following
two input parameters at runtime:

::

  integration.type = SUNDIALS
  integration.sundials.strategy = ERK

If instead one wishes to use the SUNDIALS multirate integrator, then the user
will need to use the following runtime inputs parameters:

::

  integration.type = SUNDIALS
  integration.sundials.strategy = MRI

In addition, to set up the multirate problem, the user needs to supply a fast
timescale right-hand-side function in addition to the usual right hand side
function (which is interpreted as the slow timescale right-hand side). The user
will also need to supply the ratio of the slow timestep size to the fast
timestep size, which is an integer corresponding to the number of fast
timesteps the integrator will take per every slow timestep. An example code
snippet would look as follows:

.. highlight:: c++

::

   #include <AMReX_TimeIntegrator.H>

   Vector<MultiFab> Sborder; // Vector of MultiFab(s) containing old-time state data and ghost cells
   Vector<MultiFab> Snew;    // Vector of MultiFab(s) where we want new-time state data
   Geometry geom;    // The domain (or level) geometry

   // [Fill Sborder here]

   // Create a time integrator that will work with
   // MultiFabs with the same BoxArray, DistributionMapping,
   // and number of components as the state_data MultiFab.
   TimeIntegrator<Vector<MultiFab> > integrator(Sborder);

   // Create a slow timescale RHS function we will integrate
   auto rhs_fun = [&](Vector<MultiFab>& rhs, const Vector<MultiFab>& state, const Real time){
       // User function to calculate the rhs MultiFab given the state MultiFab(s)
       fill_rhs(rhs, state, time);
   };

   // Create a fast timescale RHS function to integrate
   auto rhs_fun_fast = [&](Vector<MultiFab>& rhs,
                           const Vector<MultiFab>& stage_data,
                           const Vector<MultiFab>& state, const Real time) {
        // User function to calculate the fast-timescale rhs MultiFab given
        // the state MultiFab and stage_data which holds the previously
        // accessed slow-timescale stage state data.
        fill_fast_rhs(rhs, stage_data, state, time);
   };

   // The post update function is called after updating state data or
   // immediately before using state data to calculate a fast or slow right hand side.
   // (it is a good place to e.g. fill boundary conditions)
   auto post_update_fun = [&](Vector<MultiFab>& S_data, const Real time) {
       // Call user function to update state MultiFab(s), e.g. fill BCs
       post_update(S_data, time, geom);
   };

   // Attach the slow and fast right hand side functions to integrator
   integrator.set_rhs(rhs_fun);
   integrator.set_fast_rhs(rhs_fun_fast);

   // This sets the ratio of slow timestep size to fast timestep size as an integer,
   // or equivalently, the number of fast timesteps per slow timestep.
   integrator.set_slow_fast_timestep_ratio(2);

   // Attach the post update function to the integrator
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

The full set of integrator options are detailed as follows:

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

  ## *** Parameters Specific to SUNDIALS ERK Strategy ***
  ## (Requires integration.type=SUNDIALS and integration.sundials.strategy=ERK)
  ## integration.sundials.erk.method specifies which explicit Runge Kutta method
  ## for SUNDIALS to use. The following options are supported:
  ## "SSPRK3" = 3rd order strong stability preserving RK (default)
  ## "Trapezoid" = 2nd order trapezoidal rule
  ## "ForwardEuler" = 1st order forward euler
  ## for example:
  integration.sundials.erk.method = SSPRK3

  ## *** Parameters Specific to SUNDIALS MRI Strategy ***
  ## (Requires integration.type=SUNDIALS and integration.sundials.strategy=MRI)
  ## integration.sundials.mri.implicit_inner specifies whether or not to use an implicit inner solve
  ## integration.sundials.mri.outer_method specifies which outer (slow) method to use
  ## integration.sundials.mri.inner_method specifies which inner (fast) method to use
  ## The following options are supported for both the inner and outer methods:
  ## "KnothWolke3" = 3rd order Knoth-Wolke method (default for outer method)
  ## "Trapezoid" = 2nd order trapezoidal rule
  ## "ForwardEuler" = 1st order forward euler (default for inner method)
  ## for example:
  integration.sundials.mri.implicit_inner = false
  integration.sundials.mri.outer_method = KnothWolke3
  integration.sundials.mri.inner_method = Trapezoid
