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

        AMREX_D_TERM(int n_cell_x = 128;,
                     int n_cell_y = 32;,
                     int n_cell_z = 64);

        AMREX_D_TERM(int max_grid_size_x = 32;,
                     int max_grid_size_y = 32;,
                     int max_grid_size_z = 32);

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

        MultiFab mf(ba,dm,1,0);
        auto const& ma = mf.arrays();
        ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
            ma[b](i,j,k) = std::exp(-10._rt*
                (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
        });

        MultiFab mf2(ba,dm,1,0);

        // For each dimension, there are 5 possibilities
        constexpr int ncases = 5;
        Array<std::pair<FFT::Boundary,FFT::Boundary>,ncases>
            bcs{std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::periodic,
                                                       FFT::Boundary::periodic},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::even,
                                                       FFT::Boundary::even},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::even,
                                                       FFT::Boundary::odd },
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::odd,
                                                       FFT::Boundary::even},
                std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::odd,
                                                       FFT::Boundary::odd }};
        int ncasesy = (AMREX_SPACEDIM > 1) ? ncases : 1;
        int ncasesz = (AMREX_SPACEDIM > 2) ? ncases : 1;
        for (int zcase = 0; zcase < ncasesz; ++zcase) {
        for (int ycase = 0; ycase < ncasesy; ++ycase) {
        for (int xcase = 0; xcase < ncases ; ++xcase) {
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
                fft_bc{AMREX_D_DECL(bcs[xcase],bcs[ycase],bcs[zcase])};
            amrex::Print() << "  Testing (";
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                amrex::Print() << "(" << getEnumNameString(fft_bc[idim].first)
                               << "," << getEnumNameString(fft_bc[idim].second)
                               << ")";
                if (idim+1 < AMREX_SPACEDIM) { amrex::Print() << " "; }
            }
            amrex::Print() << ")\n";

            mf2.setVal(std::numeric_limits<Real>::max());

            FFT::R2X fft(geom.Domain(), fft_bc);
            auto scaling_factor = fft.scalingFactor();
            fft.forwardThenBackward(mf, mf2, [=] AMREX_GPU_DEVICE (int, int, int, auto& sp)
            {
                sp *= scaling_factor;
            });

            MultiFab::Subtract(mf2, mf, 0, 0, 1, 0);

            auto error = mf2.norminf();
            amrex::Print() << "    Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-6f;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }}}
    }
    amrex::Finalize();
}
