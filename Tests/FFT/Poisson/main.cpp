#include <AMReX_FFT_Poisson.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void make_rhs (MultiFab& rhs, Geometry const& geom,
               Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> const& fft_bc)
{
    auto const& dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
    auto const& probhi = geom.ProbHiArray();
    GpuArray<Real,AMREX_SPACEDIM> center
        {AMREX_D_DECL(0.5_rt*(problo[0]+probhi[0]),
                      0.5_rt*(problo[1]+probhi[1]),
                      0.5_rt*(problo[2]+probhi[2]))};
    GpuArray<Real,AMREX_SPACEDIM> problen
        {AMREX_D_DECL((probhi[0]-problo[0]),
                      (probhi[1]-problo[1]),
                      (probhi[2]-problo[2]))};

    GpuArray<Real,AMREX_SPACEDIM> fac
        {AMREX_D_DECL(2._rt*Math::pi<Real>()/problen[0],
                      2._rt*Math::pi<Real>()/problen[1],
                      2._rt*Math::pi<Real>()/problen[2])};

    auto const& rhsma = rhs.arrays();
    ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect iv(AMREX_D_DECL(i,j,k));
        Real r = 1.0_rt;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Real x = (iv[idim]+0.5_rt) * dx[idim];
            if (fft_bc[idim].first == FFT::Boundary::periodic) {
                r *= (0.11_rt + std::sin((x+0.1_rt)*fac[idim]));
            } else if (fft_bc[idim].first == FFT::Boundary::even &&
                       fft_bc[idim].second == FFT::Boundary::even) {
                r *= (0.12_rt + std::cos(x*2._rt*fac[idim]));
            } else if (fft_bc[idim].first == FFT::Boundary::odd &&
                       fft_bc[idim].second == FFT::Boundary::odd) {
                r *= std::sin(x*1.5_rt*fac[idim]);
            } else if (fft_bc[idim].first == FFT::Boundary::odd &&
                               fft_bc[idim].second == FFT::Boundary::even) {
                r *= std::sin(x*0.75_rt*fac[idim]);
            } else if (fft_bc[idim].first == FFT::Boundary::even &&
                       fft_bc[idim].second == FFT::Boundary::odd) {
                r *= std::cos(x*0.75_rt*fac[idim]);
            }
            x -= center[idim];
            x /= problen[idim];
            r *= 1.0_rt + 0.1_rt*Math::abs(std::tanh(x));
        }
        rhsma[b](i,j,k) = r;
    });

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
}

std::pair<Real,Real> check_convergence
    (MultiFab const& phi, MultiFab const& rhs, Geometry const& geom)
{
    MultiFab res(phi.boxArray(), phi.DistributionMap(), 1, 0);
    auto const& res_ma = res.arrays();
    auto const& phi_ma = phi.const_arrays();
    auto const& rhs_ma = rhs.const_arrays();
    auto const& dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> lapfac
        {AMREX_D_DECL(1._rt/(dx[0]*dx[0]),
                      1._rt/(dx[1]*dx[1]),
                      1._rt/(dx[2]*dx[2]))};
    ParallelFor(res, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        auto const& phia = phi_ma[b];
        Real lap = (phia(i-1,j,k)-2._rt*phia(i,j,k)+phia(i+1,j,k)) * lapfac[0];
#if (AMREX_SPACEDIM >= 2)
        lap += (phia(i,j-1,k)-2._rt*phia(i,j,k)+phia(i,j+1,k)) * lapfac[1];
#endif
#if (AMREX_SPACEDIM == 3)
        lap += (phia(i,j,k-1)-2._rt*phia(i,j,k)+phia(i,j,k+1)) * lapfac[2];
#endif
        res_ma[b](i,j,k) = rhs_ma[b](i,j,k) - lap;
    });
    auto bnorm = rhs.norminf();
    auto rnorm = res.norminf();
    return {bnorm, rnorm};
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        AMREX_D_TERM(int n_cell_x = 64;,
                     int n_cell_y = 48;,
                     int n_cell_z = 128);

        AMREX_D_TERM(int max_grid_size_x = 32;,
                     int max_grid_size_y = 32;,
                     int max_grid_size_z = 32);

        AMREX_D_TERM(Real prob_lo_x = 0.;,
                     Real prob_lo_y = 0.;,
                     Real prob_lo_z = 0.);
        AMREX_D_TERM(Real prob_hi_x = 1.1;,
                     Real prob_hi_y = 0.8;,
                     Real prob_hi_z = 1.9);

        {
            ParmParse pp;
            AMREX_D_TERM(pp.query("n_cell_x", n_cell_x);,
                         pp.query("n_cell_y", n_cell_y);,
                         pp.query("n_cell_z", n_cell_z));
            AMREX_D_TERM(pp.query("max_grid_size_x", max_grid_size_x);,
                         pp.query("max_grid_size_y", max_grid_size_y);,
                         pp.query("max_grid_size_z", max_grid_size_z));
        }

        Box domain(IntVect(0),IntVect(AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)));
        BoxArray ba(domain);
        ba.maxSize(IntVect(AMREX_D_DECL(max_grid_size_x,
                                        max_grid_size_y,
                                        max_grid_size_z)));
        DistributionMapping dm(ba);

        Geometry geom;
        {
            geom.define(domain,
                        RealBox(AMREX_D_DECL(prob_lo_x,prob_lo_y,prob_lo_z),
                                AMREX_D_DECL(prob_hi_x,prob_hi_y,prob_hi_z)),
                        CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        }

        // For each dimension, there are 5 possibilities
        constexpr int ncases = 5;
        Array<std::pair<FFT::Boundary,FFT::Boundary>,ncases>
            bcs{std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::periodic,
                                                       FFT::Boundary::periodic},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::odd,
                                                       FFT::Boundary::odd},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::even,
                                                       FFT::Boundary::even},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::odd,
                                                       FFT::Boundary::even},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::even,
                                                       FFT::Boundary::odd}};

        int ncasesy = (AMREX_SPACEDIM > 1) ? ncases : 1;
        int ncasesz = (AMREX_SPACEDIM > 2) ? ncases : 1;
        int icase = 0;
        for (int zcase = 0; zcase < ncasesz; ++zcase) {
        for (int ycase = 0; ycase < ncasesy; ++ycase) {
        for (int xcase = 0; xcase < ncases ; ++xcase) {
            ++icase;
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
                fft_bc{AMREX_D_DECL(bcs[xcase],bcs[ycase],bcs[zcase])};
            amrex::Print() << "  (" << icase << ") Testing (";
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                amrex::Print() << "(" << getEnumNameString(fft_bc[idim].first)
                               << "," << getEnumNameString(fft_bc[idim].second)
                               << ")";
                if (idim+1 < AMREX_SPACEDIM) { amrex::Print() << " "; }
            }
            amrex::Print() << ")\n";

            MultiFab rhs(ba,dm,1,0);
            MultiFab soln(ba,dm,1,1);
            soln.setVal(std::numeric_limits<Real>::max());
            make_rhs(rhs, geom, fft_bc);

            FFT::Poisson fft_poisson(geom, fft_bc);
            fft_poisson.solve(soln, rhs);

            auto [bnorm, rnorm] = check_convergence(soln, rhs, geom);
            amrex::Print() << "       rhs inf norm " << bnorm << "\n"
                           << "       res inf norm " << rnorm << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 2.e-3f;
#else
            auto eps = 1.e-11;
#endif
            AMREX_ALWAYS_ASSERT(rnorm < eps*bnorm);
        }}}

#if (AMREX_SPACEDIM == 3)
        amrex::Print() << "  Testing PoissonHybrid\n";

        icase = 0;
        for (int ycase = 0; ycase < ncasesy; ++ycase) {
        for (int xcase = 0; xcase < ncases ; ++xcase) {
            ++icase;
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
                fft_bc{bcs[xcase], bcs[ycase],
                       std::make_pair(FFT::Boundary::even,FFT::Boundary::even)};
            amrex::Print() << "  (" << icase << ") Testing (";
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                amrex::Print() << "(" << getEnumNameString(fft_bc[idim].first)
                               << "," << getEnumNameString(fft_bc[idim].second)
                               << ")";
                if (idim+1 < AMREX_SPACEDIM) { amrex::Print() << " "; }
            }
            amrex::Print() << ")\n";

            MultiFab rhs(ba,dm,1,0);
            MultiFab soln(ba,dm,1,1);
            soln.setVal(std::numeric_limits<Real>::max());
            make_rhs(rhs, geom, fft_bc);

            Gpu::DeviceVector<Real> dz(n_cell_z, geom.CellSize(2));
            // or Vector<Real> dz(n_cell_z, geom.CellSize(2));

            FFT::PoissonHybrid fft_poisson(geom, fft_bc);
            fft_poisson.solve(soln, rhs, dz);

            auto [bnorm, rnorm] = check_convergence(soln, rhs, geom);
            amrex::Print() << "       rhs inf norm " << bnorm << "\n"
                           << "       res inf norm " << rnorm << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 2.e-3f;
#else
            auto eps = 1.e-11;
#endif
            AMREX_ALWAYS_ASSERT(rnorm < eps*bnorm);
        }}
#endif
    }

    amrex::Finalize();
}
