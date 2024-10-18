#include <AMReX_FFT.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        AMREX_D_TERM(int n_cell_x = 64;,
                     int n_cell_y = 64;,
                     int n_cell_z = 64);

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
        }

        Box domain(IntVect(0),IntVect(AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)));
        BoxArray ba = amrex::decompose(domain, ParallelDescriptor::NProcs());
        DistributionMapping dm(ba);

        Geometry geom;
        {
            geom.define(domain,
                        RealBox(AMREX_D_DECL(prob_lo_x,prob_lo_y,prob_lo_z),
                                AMREX_D_DECL(prob_hi_x,prob_hi_y,prob_hi_z)),
                        CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        }
        auto const& dx = geom.CellSizeArray();

        MultiFab rhs(ba,dm,1,0);
        MultiFab soln(ba,dm,1,0);
        auto const& rhsma = rhs.arrays();
        ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
            rhsma[b](i,j,k) = std::exp(-10._rt*
                (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
        });

        // Shift rhs so that its sum is zero.
        auto rhosum = rhs.sum(0);
        rhs.plus(-rhosum/geom.Domain().d_numPts(), 0, 1);

        auto t0 = amrex::second();

        FFT::R2C fft(geom.Domain());
        FFT::PoissonSpectral<Real> post_forward(geom);

        auto t1 = amrex::second();

        double tsolve;
        for (int n = 0; n < 2; ++n) {
            auto ta = amrex::second();
            fft.forwardThenBackward(rhs, soln, post_forward);
            auto tb = amrex::second();
            tsolve = tb-ta;
        }

        amrex::Print() << "  AMReX FFT setup time: " << t1-t0 << ", solve time "
                       << tsolve << "\n";

        {
            MultiFab phi(soln.boxArray(), soln.DistributionMap(), 1, 1);
            MultiFab res(soln.boxArray(), soln.DistributionMap(), 1, 0);
            MultiFab::Copy(phi, soln, 0, 0, 1, 0);
            phi.FillBoundary(geom.periodicity());
            auto const& res_ma = res.arrays();
            auto const& phi_ma = phi.const_arrays();
            auto const& rhs_ma = rhs.const_arrays();
            ParallelFor(res, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& phia = phi_ma[b];
                auto lap = AMREX_D_TERM
                    (((phia(i-1,j,k)-2.*phia(i,j,k)+phia(i+1,j,k)) / (dx[0]*dx[0])),
                   + ((phia(i,j-1,k)-2.*phia(i,j,k)+phia(i,j+1,k)) / (dx[1]*dx[1])),
                   + ((phia(i,j,k-1)-2.*phia(i,j,k)+phia(i,j,k+1)) / (dx[2]*dx[2])));
                res_ma[b](i,j,k) = rhs_ma[b](i,j,k) - lap;
            });
            auto bnorm = rhs.norminf();
            auto rnorm = res.norminf();
            amrex::Print() << "  rhs inf norm " << bnorm << "\n"
                           << "  res inf norm " << rnorm << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-3f;
#else
            auto eps = 1.e-11;
#endif
            AMREX_ALWAYS_ASSERT(rnorm < eps*bnorm);
        }
    }
    amrex::Finalize();
}
