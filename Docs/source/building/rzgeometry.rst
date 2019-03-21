Building WarpX to use RZ geometry
=================================

WarpX can be built to run with RZ geometry. Currently, this only allows pure axisymmetry (i.e. mode 0) with an FDTD solver.

To select RZ geometry, set the flag USE_RZ = TRUE when compiling:
::

    make -j 4 USE_RZ=TRUE

The executable produced will have "RZ" as a suffix.
