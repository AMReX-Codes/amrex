#include <AMReX_FFT_Poisson.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

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
        auto const& dx = geom.CellSizeArray();
        GpuArray<Real,AMREX_SPACEDIM> center
            {AMREX_D_DECL(0.5_rt*(prob_lo_x+prob_hi_x),
                          0.5_rt*(prob_lo_y+prob_hi_y),
                          0.5_rt*(prob_lo_z+prob_hi_z))};
        GpuArray<Real,AMREX_SPACEDIM> problen
            {AMREX_D_DECL((prob_hi_x-prob_lo_x),
                          (prob_hi_y-prob_lo_y),
                          (prob_hi_z-prob_lo_z))};

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

            GpuArray<Real,AMREX_SPACEDIM> fac
                {AMREX_D_DECL(2._rt*Math::pi<Real>()/problen[0],
                              2._rt*Math::pi<Real>()/problen[1],
                              2._rt*Math::pi<Real>()/problen[2])};

            MultiFab rhs(ba,dm,1,0);
            MultiFab soln(ba,dm,1,0);
            soln.setVal(std::numeric_limits<Real>::max());

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

            // We know that the sum of our rhs is zero for non-Dirichlet
            // cases. Otherwise, we should shift rhs so that its sum is zero.

            FFT::Poisson fft_poisson(geom, fft_bc);
            fft_poisson.solve(soln, rhs);

            MultiFab phi(soln.boxArray(), soln.DistributionMap(), 1, 1);
            MultiFab res(soln.boxArray(), soln.DistributionMap(), 1, 0);
            MultiFab::Copy(phi, soln, 0, 0, 1, 0);
            phi.FillBoundary(geom.periodicity());
            auto const& res_ma = res.arrays();
            auto const& phi_ma = phi.const_arrays();
            auto const& rhs_ma = rhs.const_arrays();
            GpuArray<Real,AMREX_SPACEDIM> lapfac
                {AMREX_D_DECL(1._rt/(dx[0]*dx[0]),
                              1._rt/(dx[1]*dx[1]),
                              1._rt/(dx[2]*dx[2]))};
            ParallelFor(res, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& phia = phi_ma[b];
                Real lap = 0;
                if (i == 0 && fft_bc[0].first == FFT::Boundary::odd) {
                    lap += (-3._rt*phia(i,j,k)+phia(i+1,j,k)) * lapfac[0];
                } else if (i == 0 && fft_bc[0].first == FFT::Boundary::even) {
                    lap += (-phia(i,j,k)+phia(i+1,j,k)) * lapfac[0];
                } else if (i == n_cell_x-1 && fft_bc[0].second == FFT::Boundary::odd) {
                    lap += (phia(i-1,j,k)-3._rt*phia(i,j,k)) * lapfac[0];
                } else if (i == n_cell_x-1 && fft_bc[0].second == FFT::Boundary::even) {
                    lap += (phia(i-1,j,k)-phia(i,j,k)) * lapfac[0];
                } else {
                    lap += (phia(i-1,j,k)-2._rt*phia(i,j,k)+phia(i+1,j,k)) * lapfac[0];
                }
#if (AMREX_SPACEDIM >= 2)
                if (j == 0 && fft_bc[1].first == FFT::Boundary::odd) {
                    lap += (-3._rt*phia(i,j,k)+phia(i,j+1,k)) * lapfac[1];
                } else if (j == 0 && fft_bc[1].first == FFT::Boundary::even) {
                    lap += (-phia(i,j,k)+phia(i,j+1,k)) * lapfac[1];
                } else if (j == n_cell_y-1 && fft_bc[1].second == FFT::Boundary::odd) {
                    lap += (phia(i,j-1,k)-3._rt*phia(i,j,k)) * lapfac[1];
                } else if (j == n_cell_y-1 && fft_bc[1].second == FFT::Boundary::even) {
                    lap += (phia(i,j-1,k)-phia(i,j,k)) * lapfac[1];
                } else {
                    lap += (phia(i,j-1,k)-2._rt*phia(i,j,k)+phia(i,j+1,k)) * lapfac[1];
                }
#endif
#if (AMREX_SPACEDIM == 3)
                if (k == 0 && fft_bc[2].first == FFT::Boundary::odd) {
                    lap += (-3._rt*phia(i,j,k)+phia(i,j,k+1)) * lapfac[2];
                } else if (k == 0 && fft_bc[2].first == FFT::Boundary::even) {
                    lap += (-phia(i,j,k)+phia(i,j,k+1)) * lapfac[2];
                } else if (k == n_cell_z-1 && fft_bc[2].second == FFT::Boundary::odd) {
                    lap += (phia(i,j,k-1)-3._rt*phia(i,j,k)) * lapfac[2];
                } else if (k == n_cell_z-1 && fft_bc[2].second == FFT::Boundary::even) {
                    lap += (phia(i,j,k-1)-phia(i,j,k)) * lapfac[2];
                } else {
                    lap += (phia(i,j,k-1)-2._rt*phia(i,j,k)+phia(i,j,k+1)) * lapfac[2];
                }
#endif
                res_ma[b](i,j,k) = rhs_ma[b](i,j,k) - lap;
            });
            auto bnorm = rhs.norminf();
            auto rnorm = res.norminf();
            amrex::Print() << "       rhs inf norm " << bnorm << "\n"
                           << "       res inf norm " << rnorm << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 2.e-3f;
#else
            auto eps = 1.e-11;
#endif
            AMREX_ALWAYS_ASSERT(rnorm < eps*bnorm);
        }}}
    }
    amrex::Finalize();
}
