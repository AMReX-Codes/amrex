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

Please refer to the AMReX documentation at :ref:`amrex_docs:swfftdoc` for a brief explanation of how the SWFFT redistributes data into pencil grids.

AMReX contains two SWFFT tutorials, ``SWFFT_poisson`` and ``SWFFT_simple``:

- ``SWFFT_poisson`` tutorial: The tutorial found in ``amrex/Tutorials/SWFFT/SWFFT_poisson`` solves a Poisson equation with periodic boundary conditions.  In it, both a forward FFT and reverse FFT are called to solve the equation, however, no reordering of the DFT data in k-space is performed.

- ``SWFFT_simple`` tutorial: This tutorial: ``amrex/Tutorials/SWFFT/SWFFT_simple``, is useful if the objective is to simply take a forward FFT of data, and the DFT's ordering in k-space matters to the user.  This tutorial initializes a 3D or 2D :cpp:`MultiFab`, takes a forward FFT, and then redistributes the data in k-space back to the "correct," 0 to :math:`2\pi`, ordering.  The results are written to a plot file.

.. toctree::
   :maxdepth: 1

.. _section:swfft_tutorial:swfft_pois:

SWFFT_poisson
--------------------------

In this test case we set up a right hand side (rhs), call the forward transform,
modify the coefficients, then call the backward solver and output the solution
to the discrete Poisson equation.

To build the code, type 'make' in ``amrex/Tutorials/SWFFT/SWFFT_poisson``.  This
will include code from ``amrex/Src/Extern/SWFFT`` and you will need to
link to the FFT solvers themselves (on NERSC's Cori machine, for example,
you would need to "module load fft")

To run the code, type 'main3d.gnu.MPI.ex inputs' in this directory

To visualize the output, set the bool write_data to true, then
use amrvis3d (source available at https://github.com/AMReX-Codes/Amrvis):

amrvis3d -mf RHS SOL_EXACT SOL_COMP

to visualize the rhs, the exact solution and the computed solution.

The max norm of the difference between the exact and computed solution is also printed.

For instructions on how to take a forward FFT only using SWFFT, please refer to :ref:`section:swfft_tutorial:swfft_simple`.

.. _section:swfft_tutorial:swfft_simple:

SWFFT_simple
--------------------------

This tutorial initializes a 3D or 2D :cpp:`MultiFab`, takes a forward FFT, and then redistributes the data in k-space back to the "correct," 0 to :math:`2\pi`, ordering.  The results are written to a plot file.

In a similar fashion to the ``SWFFT_poisson`` tutorial:

To build the code, type 'make' in ``amrex/Tutorials/SWFFT/SWFFT_simple``.  This
will include code from ``amrex/Src/Extern/SWFFT`` and you will need to
link to the FFT solvers themselves (on NERSC's Cori machine, for example,
you would need to "module load fft")

To run the code, type 'main*.ex inputs.oneGrid' in this directory to run the code in serial. To run the code in parallel, type 'mpiexec -n $N main*.ex inputs.multipleGrids' instead, where ``N`` holds the number of MPI processes (equal to the number of grids). ``run_me_2d`` and ``run_me_3d`` also provide examples of how to run the code.

Use amrvis2d or amrvis3d to visualize the output (source available at https://github.com/AMReX-Codes/Amrvis):

amrvis${dims}d plt_fft*

where ``dims`` specifies ``AMREX_SPACEDIM``. The DFT of the data and the original data are labeled as ``FFT_of_phi`` and ``phi`` within the plot file.

The :ref:`section:swfft_tutorial:swfft_pois` tutorial provides an example of solving a Poisson equation using a discrete spectral method, in which a forward and reverse FFT of a :cpp:`MultiFab` are computed.
