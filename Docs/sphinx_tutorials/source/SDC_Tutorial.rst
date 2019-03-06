.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/SDC
==========================

**MISDC_ADR_2d**
----------------

This tutorial presents an example of using a "multi-implicit" spectral
deferred corrections (MISDC) integrator to solve a simple scalar
advection-diffusion-reaction equation in two dimensions.  Both
diffusion and reaction terms are treated implicitly but solved for
independently in an operator splitting fashion.  The advection is
treated explicitly.  The relative strengths of the three terms can be
adjusted by changing the coefficients a,d, and r in inputs_2d.

The advection operator is a 4th-order centered difference in flux
form.  The diffusion operator is a 2nd order discretization of the
Laplacian, and the implicit diffusion solve is done using
multigrid. The "reaction" term here is just a simple linear damping
hence the implicit solve is trivial.  See the routines in
functions_2d.f90 for the code that evaluates the rhs terms.

The simple form of the equation allows for an exact solution of the
PDE in a periodic geometry. There is a flag called "plot_err" in
main.cpp, which if set equal 1 will cause the code to output the error
in the solution for plotting.  If the advection term is omitted (a=0),
then an exact solution to the method of lines ODE is computed and used
to compute the error.  Hence the error in this case will scale in dt
with the order of the time integrator.

This code can also be run as an IMEX advection-diffusion example
simply by setting Npieces=2 in main.cpp.  This should also be
equivalent to setting r=0.
