
To build `libwarpx.so` for a test problem say Langmuir:

  `cd ../Langmuir`
  `make -j4 USE_PYTHON_MAIN=TRUE`

To use `libwarpx.so`, copy it to this directory and do

  `make`

We now have an executable.  To run the test,

  `mpiexec -n 4 ./main3d.Linux.gcc.gfortran.MPI.ex ../Langmuir/inputs`
