Building WarpX to use RZ geometry
=================================

WarpX can be built to run with RZ geometry. Both an FDTD solver (the default)
and a PSATD solver are available. Both solvers allow multiple azimuthal modes.

To select RZ geometry, set the flag USE_RZ = TRUE when compiling:
::

    make -j 4 USE_RZ=TRUE

Note that this sets DIM=2, which is required with USE_RZ=TRUE.
The executable produced will have "RZ" as a suffix.

RZ geometry with spectral solver
--------------------------------

Additional steps are needed to build the spectral solver. Some of the steps
are the same as is done for the Cartesian spectral solver, setting up the FFTW
package and setting ``USE_PSATD=TRUE``.

      - Install (or load) an MPI-enabled version of FFTW.
        For instance, for Debian, this can be done with
        ::

           apt-get install libfftw3-dev libfftw3-mpi-dev

      - Set the environment variable ``FFTW_HOME`` to the path for FFTW.
        For instance, for Debian, this is done with
        ::

           export FFTW_HOME=/usr/

      - Download and build the blaspp and lapackpp packages. These can be obtained from bitbucket using mercurial.
        ::

           hg clone https://bitbucket.org/icl/blaspp
           hg clone https://bitbucket.org/icl/lapackpp

      - Set the environment variables BLASPP_HOME and LAPACKPP_HOME to the locations where
        the packages were downloaded and built. For example, using bash:
        ::

           export BLASPP_HOME=/location/of/packages/blaspp
           export LAPACKPP_HOME=/location/of/packages/lapackpp

      - Set ``USE_PSATD=TRUE`` when compiling:
        ::

           make -j 4 USE_RZ=TRUE USE_PSATD=TRUE
