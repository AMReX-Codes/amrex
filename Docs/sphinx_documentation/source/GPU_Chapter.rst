.. _Chap:GPU:

GPU
===

In this chapter, we will present the GPU support in AMReX.  AMReX targets
NVIDIA, AMD and Intel GPUs using their native vendor language and therefore
requires CUDA, HIP/ROCm and DPC++/SYCL, for NVIDIA, AMD and Intel GPUs, respectively.
Users can also use OpenMP and/or OpenACC in their applications.

AMReX supports NVIDIA GPUs with compute capability >= 6 and CUDA >= 10
as well as CUDA 9.1.  While HIP and DPC++ compilers are in development in
preparation for Frontier and Aurora, AMReX only supports the latest
publicly released versions of those compilers on the Iris and Tulip testbeds.

For complete details of CUDA, HIP, DPC++, OpenMP and OpenACC
languages, see their respective documentations.

Be aware, this documentation is currently focused on CUDA.  HIP and DPC++ documentation
is forthcoming.

A number of tutorials can be found at `Tutorials/GPU`_.

.. _`Tutorials/GPU`: https://amrex-codes.github.io/amrex/tutorials_html/GPU_Tutorial.html

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   GPU
