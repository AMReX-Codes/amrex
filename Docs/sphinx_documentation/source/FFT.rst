.. role:: cpp(code)
   :language: c++

.. _sec:FFT:r2c:

FFT::R2C Class
==============

Class template `FFT::R2C` supports discrete Fourier transforms between real
and complex data. The name R2C indicates that the forward transform converts
real data to complex data, while the backward transform converts complex
data to real data. It should be noted that both directions of transformation
are supported, not just from real to complex.

The implementation utilizes cuFFT, rocFFT, oneMKL and FFTW, for CUDA, HIP,
SYCL and CPU builds, respectively. Because the parallel communication is
handled by AMReX, it does not need the parallel version of
FFTW. Furthermore, there is no constraint on the domain decomposition such
as one Box per process. This class performs parallel FFT on AMReX's parallel
data containers (e.g., :cpp:`MultiFab` and
:cpp:`FabArray<BaseFab<ComplexData<Real>>>`. For local FFT, the users can
use FFTW, cuFFT, rocFFT, or oneMKL directly.

Other than using column-majored order, AMReX follows the convention of
FFTW. Applying the forward transform followed by the backward transform
scales the original data by the size of the input array. The layout of the
complex data also follows the FFTW convention, where the complex Hermitian
output array has `(nx/2+1,ny,nz)` elements. Here `nx`, `ny` and `nz` are the
sizes of the real array and the division is rounded down.

Below are examples of using :cpp:`FFT:R2C`.

.. highlight:: c++

::

    Geometry geom(...);
    MultiFab mfin(...);
    MultiFab mfout(...);

    auto scaling = 1. / geom.Domain().d_numPts();

    FFT::R2C r2c(geom.Domain());
    r2c.forwardThenBackward(mfin, mfout,
        [=] AMREX_GPU_DEVICE (int, int, int, auto& sp)
        {
            sp *= scaling;
        });

    cMultiFab cmf(...);
    FFT::R2C<Real,FFT::Direction::forward> r2c_forward(geom.Domain());
    r2c_forward(mfin, cmf);

    FFT::R2C<Real,FFT::Direction::backward> r2c_backward(geom.Domain());
    r2c_backward(cmf, mfout);

Note that using :cpp:`forwardThenBackward` is expected to be more efficient
than separate calls to :cpp:`forward` and :cpp:`backward` because some
parallel communication can be avoided. It should also be noted that a lot of
preparation works are done in the construction of an :cpp:`FFT::R2C`
object. Therefore, one should cache it for reuse if possible.


Poisson Solver
==============

AMReX provides FFT based Poisson solvers. :cpp:`FFT::Poisson` supports all
periodic boundaries using purely FFT. :cpp:`FFT::PoissonHybrid` is a 3D only
solver that supports periodic boundaries in the first two dimensions and
Neumann boundary in the last dimension. Similar to :cpp:`FFT::R2C`, the
Poisson solvers should be cached for reuse.
