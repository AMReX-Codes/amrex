
For CUDA C++ only, you can compile this tutorial with

    make COMP=gcc
or other compilers.

For CUDA C++ and OpenACC, compile with

    make COMP=pgi USE_ACC=TRUE

For CUDA C++ and OpenMP, compile with

    make COMP=pgi USE_OMP_OFFLOAD=TRUE
