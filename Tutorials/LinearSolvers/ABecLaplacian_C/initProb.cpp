
#include "MyTest.H"
#include "initProb_K.H"

using namespace amrex;

void
MyTest::initProbPoisson ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        const auto prob_lo = geom[ilev].ProbLoArray();
        const auto dx      = geom[ilev].CellSizeArray();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto rhsfab = rhs[ilev].array(mfi);
            auto exactfab = exact_solution[ilev].array(mfi);
            AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                actual_init_poisson(i,j,k,rhsfab,exactfab,prob_lo,dx);
            });
        }

        solution[ilev].setVal(0.0);
    }
}

void
MyTest::initProbABecLaplacian ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        const auto prob_lo = geom[ilev].ProbLoArray();
        const auto prob_hi = geom[ilev].ProbHiArray();
        const auto dx      = geom[ilev].CellSizeArray();
        auto a = ascalar;
        auto b = bscalar;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = mfi.growntilebox(1);

            auto rhsfab = rhs[ilev].array(mfi);
            auto exactfab = exact_solution[ilev].array(mfi);
            auto acoeffab = acoef[ilev].array(mfi);
            auto bcoeffab = bcoef[ilev].array(mfi);

            AMREX_PARALLEL_FOR_3D (gbx, i, j, k,
            {
                actual_init_bcoef(i,j,k,bcoeffab,prob_lo,prob_hi,dx);
            });

            AMREX_PARALLEL_FOR_3D (bx, i, j, k,
            {
                actual_init_abeclap(i,j,k,rhsfab,exactfab,acoeffab,bcoeffab,
                                    a,b,prob_lo,prob_hi,dx);
            });
        }

        solution[ilev].setVal(0.0);
    }
}
