.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/SWFFT
==========================

This Tutorial demonstrates how to call the SWFFT wrapper to the FFTW3 solver.

Note that the SWFFT source code was developed by Adrian Pope and colleagues and
is available at:

https://xgitlab.cels.anl.gov/hacc/SWFFT

In this test case we set up a right hand side (rhs), call the forward transform,
modify the coefficients, then call the backward solver and output the solution
to the discrete Poisson equation.

To build the code, type 'make' in amrex/Tutorials/SWFFT.  This
will include code from amrex/Src/Extern/SWFFT and you will need to
link to the FFT solvers themselves (on NERSC's Cori machine, for example,
you would need to "module load fft")

To run the code, type 'main3d.gnu.MPI.ex inputs' in this directory

To visualize the output, set the bool write_data to true, then
use amrvis3d (source available at https://github.com/AMReX-Codes/Amrvis):

amrvis3d -mf RHS SOL_EXACT SOL_COMP

to visualize the rhs, the exact solution and the computed solution.

The max norm of the difference between the exact and computed solution is also printed.


