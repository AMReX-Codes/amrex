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

Below are examples of using :cpp:`FFT::R2C`.

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

    // Use R2C provided spectral data layout.
    auto const& [cba, cdm] = r2c.getSpectralDataLayout();
    cMultiFab cmf(cba, cdm, 1, 0);
    FFT::R2C<Real,FFT::Direction::forward> r2c_forward(geom.Domain());
    r2c_forward(mfin, cmf);

    FFT::R2C<Real,FFT::Direction::backward> r2c_backward(geom.Domain());
    r2c_backward(cmf, mfout);

Note that using :cpp:`forwardThenBackward` is expected to be more efficient
than separate calls to :cpp:`forward` and :cpp:`backward` because some
parallel communication can be avoided. For the spectral data, the example
above builds :cpp:`cMultiFab` using :cpp:`FFT::R2C` provided layout. You can
also use your own :cpp:`BoxArray` and :cpp:`DistributionMapping`, but it
might result in extra communication. It should also be noted that a lot of
preparation works are done in the construction of an :cpp:`FFT::R2C`
object. Therefore, one should cache it for reuse if possible. Although
:cpp:`FFT::R2C` does not have a default constructor, one could always use
:cpp:`std::unique_ptr<FFT::R2C<Real>>` to store an object in one's class.


Poisson Solver
==============

AMReX provides FFT based Poisson solvers. Here, Poisson's equation is

.. math::

  \nabla^2 \phi = \rho.

:cpp:`FFT::Poisson` supports periodic (:cpp:`FFT::Boundary::periodic`),
homogeneous Neumann (:cpp:`FFT::Boundary::even`), and homogeneous Dirichlet
(:cpp:`FFT::Boundary::odd`) boundaries using FFT. Below is an example of
using the solver.

.. highlight:: c++

::

    Geometry geom(...);
    MultiFab soln(...);
    MultiFab rhs(...);

    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
            fft_bc{...};

    bool has_dirichlet = false;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        has_dirichlet = has_dirichlet ||
            fft_bc[idim].first == FFT::Boundary::odd ||
            fft_bc[idim].second == FFT::Boundary::odd;
    }
    if (! has_dirichlet) {
        // Shift rhs so that its sum is zero.
        auto rhosum = rhs.sum(0);
        rhs.plus(-rhosum/geom.Domain().d_numPts(), 0, 1);
    }

    FFT::Poisson fft_poisson(geom, fft_bc);
    fft_poisson.solve(soln, rhs);

:cpp:`FFT::PoissonOpenBC` is a 3D only solver that supports open
boundaries. Its implementation utilizes :cpp:`FFT::OpenBCSolver`, which can
be used for implementing convolution based solvers with a user provided
Green's function. If users want to extend the open BC solver to 2D or other
types of Green's function, they could use :cpp:`FFT::PoissonOpenBC` as an
example. Below is an example of solving Poisson's equation with open
boundaries.

.. highlight:: c++

::

    Geometry geom(...);
    MultiFab soln(...); // soln can be either nodal or cell-centered.
    MultiFab rhs(...);  // rhs must have the same index type as soln.

    int ng = ...; // ng can be non-zero, if we want to compute potential
                  // outside the domain.
    FFT::PoissonOpenBC openbc_solver(geom, soln.ixType(), IntVect(ng));
    openbc_solver.solve(soln, rhs);

:cpp:`FFT::PoissonHybrid` is a 3D only solver that supports periodic
boundaries in the first two dimensions and Neumann boundary in the last
dimension. The last dimension is solved with a tridiagonal solver that can
support non-uniform cell size in the z-direction. For most applications,
:cpp:`FFT::Poisson` should be used.

Similar to :cpp:`FFT::R2C`, the Poisson solvers should be cached for reuse,
and one might need to use :cpp:`std::unique_ptr<FFT::Poisson<MultiFab>>`
because there is no default constructor.
