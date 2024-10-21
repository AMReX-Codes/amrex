.. _Chap:FFT:

.. _sec:FFT:FFTOverview:

Discrete Fourier Transform
==========================

AMReX provides support for parallel discrete Fourier transform. The
implementation utilizes cuFFT, rocFFT, oneMKL and FFTW, for CUDA, HIP, SYCL
and CPU builds, respectively. It also provides FFT based Poisson
solvers.

.. toctree::
   :maxdepth: 1

   FFT
