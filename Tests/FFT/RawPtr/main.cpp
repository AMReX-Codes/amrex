#include <AMReX_FFT.H>

int main (int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
    amrex::Init_FFT(MPI_COMM_WORLD);
#else
    amrex::ignore_unused(argc,argv);
    amrex::Init_FFT();
#endif

    int nprocs, myproc;
#ifdef AMREX_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
    nprocs = 1;
    myproc = 0;
#endif

    using RT = amrex::Real;
    using CT = amrex::GpuComplex<RT>;

    std::array<int,AMREX_SPACEDIM> domain_size{AMREX_D_DECL(128,128,128)};

    {
        amrex::FFT::R2C<RT,amrex::FFT::Direction::both> r2c(domain_size);

        int nx = (domain_size[0] + nprocs - 1) / nprocs;
        int xlo = nx * myproc;
        nx = std::max(std::min(nx,domain_size[0]-xlo), 0);
        std::array<int,AMREX_SPACEDIM> local_start{AMREX_D_DECL(xlo,0,0)};
        std::array<int,AMREX_SPACEDIM> local_size{AMREX_D_DECL(nx,domain_size[1],domain_size[2])};

        r2c.setLocalDomain(local_start,local_size);

        auto const& [local_start_sp, local_size_sp] = r2c.getLocalSpectralDomain();
        amrex::ignore_unused(local_start_sp);

        auto nr = AMREX_D_TERM(std::size_t(local_size[0]),
                              *std::size_t(local_size[1]),
                              *std::size_t(local_size[2]));
        auto* pr = (RT*)amrex::The_Arena()->alloc(sizeof(RT)*nr);
        amrex::ParallelFor(nr, [=] AMREX_GPU_DEVICE (std::size_t i)
        {
            pr[i] = std::sin(RT(i));
        });

        auto nc = AMREX_D_TERM(std::size_t(local_size_sp[0]),
                              *std::size_t(local_size_sp[1]),
                              *std::size_t(local_size_sp[2]));
        auto* pc = (CT*)amrex::The_Arena()->alloc(sizeof(CT)*nc);

        r2c.forward(pr, pc);

        auto scaling = r2c.scalingFactor();
        amrex::ParallelFor(nc, [=] AMREX_GPU_DEVICE (std::size_t i)
        {
            pc[i] *= scaling;
        });

        r2c.backward(pc, pr);

        auto error = amrex::Reduce::Max<RT>
            (nr, [=] AMREX_GPU_DEVICE (std::size_t i)
                     {
                         return std::abs(pr[i]-std::sin(RT(i)));
                     });
        amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
        auto eps = 3.e-6F;
#else
        auto eps = 1.e-13;
#endif
        amrex::ParallelDescriptor::Barrier();
        AMREX_ALWAYS_ASSERT(error < eps);
    }

    int nbatch = 3;
    {
        amrex::FFT::Info info{};
        info.setBatchSize(nbatch);
        amrex::FFT::C2C<RT,amrex::FFT::Direction::both> c2c(domain_size,info);

        auto const& [local_start, local_size] = c2c.getLocalDomain();

        int nx = (domain_size[0] + nprocs - 1) / nprocs;
        int xlo = nx * myproc;
        nx = std::max(std::min(nx,domain_size[0]-xlo), 0);
        std::array<int,AMREX_SPACEDIM> local_start_sp{AMREX_D_DECL(xlo,0,0)};
        std::array<int,AMREX_SPACEDIM> local_size_sp{AMREX_D_DECL(nx,domain_size[1],domain_size[2])};

        c2c.setLocalSpectralDomain(local_start_sp, local_size_sp);

        auto nf = AMREX_D_TERM(std::size_t(local_size[0]),
                              *std::size_t(local_size[1]),
                              *std::size_t(local_size[2]));
        auto* pf = (CT*)amrex::The_Arena()->alloc(sizeof(CT)*nf*nbatch);
        amrex::ParallelFor(nf, [=] AMREX_GPU_DEVICE (std::size_t i)
        {
            pf[i] = amrex::GpuComplex(std::sin(RT(i)), std::cos(RT(i)));
        });

        auto nb = AMREX_D_TERM(std::size_t(local_size_sp[0]),
                              *std::size_t(local_size_sp[1]),
                              *std::size_t(local_size_sp[2]));
        auto* pb = (CT*)amrex::The_Arena()->alloc(sizeof(CT)*nb*nbatch);

        c2c.forward(pf, pb);

        auto scaling = c2c.scalingFactor();
        amrex::ParallelFor(nb, [=] AMREX_GPU_DEVICE (std::size_t i)
        {
            pb[i] *= scaling;
        });

        c2c.backward(pb, pf);

        auto error = amrex::Reduce::Max<RT>
            (nf, [=] AMREX_GPU_DEVICE (std::size_t i)
                     {
                         return amrex::norm(pf[i]-amrex::GpuComplex(std::sin(RT(i)),std::cos(RT(i))));
                     });
        amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
        auto eps = 3.e-6F;
#else
        auto eps = 1.e-13;
#endif
        amrex::ParallelDescriptor::Barrier();
        AMREX_ALWAYS_ASSERT(error < eps);
    }

    amrex::Finalize_FFT();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
