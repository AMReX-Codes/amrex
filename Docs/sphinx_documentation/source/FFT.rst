.. role:: cpp(code)
   :language: c++

.. _sec:FFT:r2c:

FFT::R2C Class
==============

Class template `FFT::R2C` supports discrete Fourier transforms between real
and complex data across MPI processes. The name R2C indicates that the
forward transform converts real data to complex data, while the backward
transform converts complex data to real data. It should be noted that both
directions of transformation are supported, not just from real to complex.

The implementation utilizes cuFFT, rocFFT, oneMKL and FFTW, for CUDA, HIP,
SYCL and CPU builds, respectively. Because the parallel communication is
handled by AMReX, it does not need the parallel version of
FFTW. Furthermore, there is no constraint on the domain decomposition such
as one Box per process. This class performs parallel FFT on AMReX's parallel
data containers (e.g., :cpp:`MultiFab` and
:cpp:`FabArray<BaseFab<ComplexData<Real>>>`.

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


Class template `FFT::R2C` also supports batched FFTs. The batch size is set
in an :cpp:`FFT::Info` object passed to the constructor of
:cpp:`FFT::R2C`. Below is an example.

.. highlight:: c++

::

    int batch_size = 10;
    Geometry geom(...);
    MultiFab mf(ba, dm, batch_size, 0);

    FFT::Info info{};
    info.setBatchSize(batch_size));
    FFT::R2C<Real,FFT::Direction::both> r2c(geom.Domain(), info);

    auto const& [cba, cdm] = r2c.getSpectralDataLayout();
    cMultiFab cmf(cba, cdm, batch_size, 0);

    r2c.forward(mf, cmf);

    // Do work on cmf.
    // Function forwardThenBackward is not yet supported for a batched FFT.

    r2c.backward(cmf, mf);

.. _sec:FFT:c2c:

FFT::C2C Class
==============

:cpp:`FFT::C2C` is a class template that supports complex to complex Fourier
transforms. It has a similar interface as :cpp:`FFT::R2C`.

.. _sec:FFT:localr2c:

FFT::LocalR2C Class
===================

Class template `FFT::LocalR2C` supports local discrete Fourier transforms
between real and complex data. The name R2C indicates that the forward
transform converts real data to complex data, while the backward transform
converts complex data to real data. It should be noted that both directions
of transformation are supported, not just from real to complex.

Below is an example of using :cpp:`FFT::LocalR2C`.

.. highlight:: c++

::

    MultiFab mf(...);
    BaseFab<GpuComplex<T>> cfab;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        FFT::LocalR2C fft(mfi.fabbox().length());
        cfab.resize(IntVect(0), fft.spectralSize()-1);
        fft.forward(mf[mfi].dataPtr(), cfab.dataPtr());
    }


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

.. _sec:FFT:rawptr:

Raw Pointer Interface
=====================

If you only want to use AMReX as a Parallel FFT library without using other
functionalities and data containers, you could use the raw pointer
interface. Below is an example.

.. highlight:: c++

::

    MPI_Init(&argc, &argv);

    // We don't need to call the full-blown amrex::Initialize
    amrex::Init_FFT(MPI_COMM_WORLD);

    int nprocs, myproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

    using RT = double;
    using CT = std::complex<RT>; // or cufftDoubleComplex, etc.

    std::array<int,3> domain_size{128,128,128};

    // FFT between real and complex.
    // Domain decomposition is flexible. The only constraint for the raw
    // pointer interface is that there can be only zero or one local box
    // per process, whereas the MultiFab interface can take any number of
    // boxes. In this case, we choose to do manual domain decomposition for
    // the real (i.e., forward) domain, and use the domain decomposition
    // provided by amrex for the complex (i.e., backward) domain.
    {
        amrex::FFT::R2C<RT,amrex::FFT::Direction::both> r2c(domain_size);

        // amrex supports 1d, 2d and 3d domain decomposition. For simplicity,
        // let's do 1d decomposition in the z-direction.
        int nz = (domain_size[2] + nprocs - 1) / nprocs;
        int zlo = nz * myproc;
        nz = std::max(std::min(nz,domain_size[2]-zlo), 0);
        std::array<int,3> local_start{0,0,zlo};
        std::array<int,3> local_size{domain_size[0],domain_size[1],nz};

        // Let amrex know the domain decomposition in the forward domain.
        r2c.setLocalDomain(local_start,local_size);

        // Use amrex's domain decomposition in the backward domain.
        auto const& [local_start_sp, local_size_sp] = r2c.getLocalSpectralDomain();

        auto nr = std::size_t(local_size[0])
                * std::size_t(local_size[1])
                * std::size_t(local_size[2]);
        auto* pr = (RT*)std::malloc(sizeof(RT)*nr); // or use cudaMalloc
        // Initialize data ...

        auto nc = std::size_t(local_size_sp[0])
                * std::size_t(local_size_sp[1])
                * std::size_t(local_size_sp[2]);
        auto* pc = (CT*)std::malloc(sizeof(CT)*nc); // or use cudaMalloc

        r2c.forward(pr, pc); // forward transform from real to complex

        // work on the complex data pointed by pc ...

        r2c.backward(pc, pr); // backward transform from complex to real

        std::free(pr);
        std::free(pc);
    }

    // Batched FFT between complex and complex.
    // In this case, we choose to use the domain decomposition provided
    // by amrex for the forward domain, and do manual domain decomposition
    // for the backward domain.
    int nbatch = 3; // batch size
    {
        amrex::FFT::Info info{};
        info.setBatchSize(nbatch);
        amrex::FFT::C2C<RT,amrex::FFT::Direction::both> c2c(domain_size,info);

        // Use amrex's domain decomposition in the forward domain.
        auto const& [local_start, local_size] = c2c.getLocalDomain();

        int nx = (domain_size[0] + nprocs - 1) / nprocs;
        int xlo = nx * myproc;
        nx = std::max(std::min(nx,domain_size[0]-xlo), 0);
        std::array<int,3> local_start_sp{xlo,0,0};
        std::array<int,3> local_size_sp{nx,domain_size[1],domain_size[2]};

        // Let amrex know the domain decomposition in the backward domain.
        c2c.setLocalSpectralDomain(local_start_sp, local_size_sp);

        auto nf = std::size_t(local_size[0])
                * std::size_t(local_size[1])
                * std::size_t(local_size[2]);
        auto* pf = (CT*)std::malloc(sizeof(CT)*nf*nbatch); // or use cudaMalloc
        // Initialize data ...

        auto nb = std::size_t(local_size_sp[0])
                * std::size_t(local_size_sp[1])
                * std::size_t(local_size_sp[2]);
        auto* pb = (CT*)std::malloc(sizeof(CT)*nb*nbatch);

        c2c.forward(pf, pb); // forward transform

        // work on the data pointed by pb

        c2c.backward(pb, pf); // backward transform

        std::free(pf);
        std::free(pb);
    }

    amrex::Finalize_FFT();

    MPI_Finalize();
