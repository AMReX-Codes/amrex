.. _Chap:GPU:

GPU
===

In this chapter, we will present the GPU support in AMReX.  Currently
AMReX only supports Nvidia GPUs.  Internally, it uses CUDA C++, but the
users can use CUDA Fortran and/or OpenACC.  A recent version of CUDA
(e.g., >= 9) is required, and the device must have compute capability
>= 6.

For complete details of CUDA, CUDA Fortran and OpenACC languages, see their
respective documentations.

A number of tutorials can be found at ``Tutorials/GPU/``.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   GPU
