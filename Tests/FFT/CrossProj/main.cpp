#include <AMReX_FFT_CrossProj.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        AMREX_D_TERM(int n_cell_x = 128;,
                     int n_cell_y = 128;,
                     int n_cell_z = 64);

        AMREX_D_TERM(int max_grid_size_x = 16;,
                     int max_grid_size_y = 16;,
                     int max_grid_size_z = 16);

        AMREX_D_TERM(Real prob_lo_x = 0.;,
                     Real prob_lo_y = 0.;,
                     Real prob_lo_z = 0.);
        AMREX_D_TERM(Real prob_hi_x = 1.;,
                     Real prob_hi_y = 1.;,
                     Real prob_hi_z = 1.);

        {
            ParmParse pp;
            AMREX_D_TERM(pp.query("n_cell_x", n_cell_x);,
                         pp.query("n_cell_y", n_cell_y);,
                         pp.query("n_cell_z", n_cell_z));
            AMREX_D_TERM(pp.query("max_grid_size_x", max_grid_size_x);,
                         pp.query("max_grid_size_y", max_grid_size_y);,
                         pp.query("max_grid_size_z", max_grid_size_z));
        }

        auto bcs = std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::periodic,
                                                          FFT::Boundary::periodic};

        Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
                fft_bc{AMREX_D_DECL(bcs,bcs,bcs)};

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

        MultiFab vel(ba,dm,AMREX_SPACEDIM,1);
        auto const& u = vel.arrays();
        ParallelFor(vel, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            constexpr Real tpi = Real(2.)*Math::pi<Real>();
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
#if (AMREX_SPACEDIM == 2)
            u[b](i,j,k,0) = std::cos(tpi*x) + std::cos(tpi*y);
            u[b](i,j,k,1) = std::cos(tpi*x+0.2_rt) + std::cos(tpi*y+0.7_rt);
#else
            u[b](i,j,k,0) = std::cos(tpi*x) + std::cos(tpi*y) + std::cos(tpi*z+0.3_rt);
            u[b](i,j,k,1) = std::cos(tpi*x+0.2_rt) + std::cos(tpi*y+0.7_rt) + std::sin(tpi*z);
            u[b](i,j,k,2) = std::sin(tpi*x+0.4_rt) + std::cos(tpi*y+0.1_rt) + std::cos(tpi*z+0.5_rt);
#endif
        });

        vel.FillBoundary(geom.periodicity());

        MultiFab vel_return(ba,dm,AMREX_SPACEDIM,1);

        {
            FFT::CrossProj crossproj(geom, fft_bc);
            crossproj.solve(vel_return,  vel ); // also fills ghost cells

            MultiFab div_orig(ba, dm, 1, 0);
            MultiFab div_err(ba, dm, 1, 0);
            auto const& ud = vel_return.arrays();
            auto const& divinit = div_orig.arrays();
            auto const& div = div_err.arrays();

            // Node-centered cross-stencil divergence at (i+1/2,j+1/2,k+1/2)
            auto nodediv = [=] AMREX_GPU_DEVICE (Array4<Real const> const& v,
                                                 int i, int j, int k) -> Real
            {
#if (AMREX_SPACEDIM == 2)
                return (v(i+1,j+1,k,0) + v(i+1,j,k,0) - v(i,j+1,k,0) - v(i,j,k,0))
                    * (0.5_rt/dx[0])
                    +  (v(i+1,j+1,k,1) + v(i,j+1,k,1) - v(i+1,j,k,1) - v(i,j,k,1))
                    * (0.5_rt/dx[1]);
#else
                Real dvx = 0, dvy = 0, dvz = 0;
                for (int s = 0; s <= 1; ++s) {
                for (int t = 0; t <= 1; ++t) {
                    dvx += v(i+1,j+s,k+t,0) - v(i,j+s,k+t,0);
                    dvy += v(i+s,j+1,k+t,1) - v(i+s,j,k+t,1);
                    dvz += v(i+s,j+t,k+1,2) - v(i+s,j+t,k,2);
                }}
                return dvx*(0.25_rt/dx[0]) + dvy*(0.25_rt/dx[1]) + dvz*(0.25_rt/dx[2]);
#endif
            };

            ParallelFor(div_err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                divinit[b](i,j,k) = nodediv(u[b], i, j, k);
                div[b](i,j,k) = nodediv(ud[b], i, j, k);
            });

            Real div_initial = div_orig.norminf();
            Real div_error = div_err.norminf();
            amrex::Print() << "  original divergence " << div_initial << "\n";
            amrex::Print() << "  divergence expected to be close to zero: " << div_error << "\n";

            // A second projection should not change the result.
            MultiFab vel_twice(ba,dm,AMREX_SPACEDIM,0);
            crossproj.solve(vel_twice, vel_return);
            MultiFab::Subtract(vel_twice, vel_return, 0, 0, AMREX_SPACEDIM, 0);
            Real idem_error = vel_twice.norminf(0, AMREX_SPACEDIM, IntVect(0));
            Real vel_max = vel_return.norminf(0, AMREX_SPACEDIM, IntVect(0));
            amrex::Print() << "  idempotency error: " << idem_error << "\n";

            Real const tol = std::is_same_v<Real,float> ? 1.e-4_rt : 1.e-12_rt;
            AMREX_ALWAYS_ASSERT(div_error <= tol*div_initial);
            AMREX_ALWAYS_ASSERT(idem_error <= tol*vel_max);
        }
    }
    amrex::Finalize();
}
