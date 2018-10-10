# Using CVODE CUDA interface with a Castro-style driver

This example program shows how to use the CUDA interface in CVODE
4.0.0 in an AMReX application.

To build the example, you should define the environment variable
`CVODE_HOME` set to the CVODE installation path (containing the
`include` and `lib` directories).

You should also add `$CVODE_HOME/lib` to `LD_LIBRARY_PATH`.

Then just do `make`.

To run, do `main3d.gnu.ex inputs`

# Building with CUDA

To build the example with CUDA CVODE, use the PGI compiler and do
`make COMP=PGI USE_CUDA=TRUE`. I've tested this with CUDA 9.

# Testing on Groot

Groot has a Xeon 5115 CPU and a NVIDIA GTX 1060 GPU. Here is the
fcompare output for the plotfiles produced by the serial CPU CVODE and
flattened CUDA CVODE approaches for a 16**3 grid:

```
            variable name            absolute error            relative error
                                        (||A - B||)         (||A - B||/||A||)
 ----------------------------------------------------------------------
 level =  1
 Y0                                0.2760985749E-05          0.2802548919E-05
 Y1                                0.1558203119E-08          0.4601430201E-04
 Y2                                0.4444774455E-05          0.3003911008E-03
```

That's pretty reasonable considering roundoff error will be different
on the CPU and GPU, they are using different linear solvers, the
linear systems are different (serialized vs flattened), and the
relative tolerance for integration is 1.e-4.

Timing for a 32**3 grid of cells:
- 57 seconds for serial CVODE
- 0.83 seconds for flattened CUDA CVODE

That's a speedup of nearly 70, which is pretty decent for a grid this
small, as the GPU can feasibly integrate many more zones than this
simultaneously.

# Using the CVODE cuSolver linear solver interface

To use the CVODE cuSolver interface, make with the flag
`USE_CVODE_CUSOLVER=TRUE` in addition to `USE_CUDA=TRUE` and define
`CUDA_HOME` to point to the CUDA toolkit installation directory
(containing `lib64`).
