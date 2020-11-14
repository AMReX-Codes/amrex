
This tutorial can be compiled with:

    make COMP=gcc

The code is mean to be run on GPU. The source file shows an example ParallelFor
instrumented with the `KernelTimer` which allows to record the thread-wise 
summed number of cycles and can be used as a proxy for GPU compute time. The 
value stored by the `KernelTimer` can be compared across kernels to give some 
quantitative idea of the relative compute work across kernels. In this example, 
the parameter `n_cells` (set in main.cpp) can be varied, and the effect on number 
of GPU cycles can be observed; for example, increasing `n_cells` by a factor of 
4 increases the domain size by a factor of 64, and is the expected factor by 
which the measured number of cycles should increase.  Here is sample output for 
cases with `n_cells=256` and `n_cells=1024`, respectively,  when this code is 
run on NVIDIA V100 GPU:

    I measured 4227657714601 cycles over all threads.
    I measured 273340022106357 cycles over all threads.
