#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        int n_cell = 128;
        int max_grid_size = 32;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        Box domain(IntVect(0), IntVect(n_cell-1));
        Periodicity period(domain.length());
        BoxArray ba(domain);
        ba.convert(IntVect(1));
        ba.maxSize(IntVect(max_grid_size));
        DistributionMapping dm{ba};

        MultiFab mf(ba, dm, 2, 2);
        {
            BoxArray ba2 = ba;
            ba2.grow(2);
            MultiFab mf2(ba2, dm, 1, 0, MFInfo{}.SetAlloc(false));
            std::unique_ptr<MultiFab> mask = mf2.OverlapMask(period);

            auto const& a = mf.arrays();
            auto const& m = mask->const_arrays();
            ParallelFor(mf, mf.nGrowVect(),
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                a[b](i,j,k,0) = Real(3.14) / m[b](i,j,k);
                a[b](i,j,k,1) = Real(6.28) / m[b](i,j,k);
            });
        }

        mf.SumBoundary(0, 1, period, false);
        mf.SumBoundary(1, 1, period, true);

        auto c0_min = mf.min(0);
        auto c0_max = mf.max(0);
        auto c1_min = mf.min(1);
        auto c1_max = mf.max(1);

        amrex::Print() << "  min(0) = " << c0_min << ", max(0) = " << c0_max
                       << ", min(1) = " << c1_min << ", max(1) = " << c1_max << '\n';

        AMREX_ALWAYS_ASSERT(amrex::almostEqual(c0_min, Real(3.14)));
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(c0_max, Real(3.14)));
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(c1_min, Real(6.28)));
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(c1_max, Real(6.28)));
    }
    amrex::Finalize();
}
