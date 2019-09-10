
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
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_bcoef(i,j,k,bcoeffab,prob_lo,prob_hi,dx);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_abeclap(i,j,k,rhsfab,exactfab,acoeffab,bcoeffab,
                                    a,b,prob_lo,prob_hi,dx);
            });
        }

        solution[ilev].setVal(0.0);
    }
}

void
MyTest::initProbABecLaplacianInhomNeumann ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        const auto prob_lo = geom[ilev].ProbLoArray();
        const auto prob_hi = geom[ilev].ProbHiArray();
        const auto dx      = geom[ilev].CellSizeArray();
        Box const& domain = geom[ilev].Domain();
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
            auto solnfab = solution[ilev].array(mfi);
            auto exactfab = exact_solution[ilev].array(mfi);
            auto acoeffab = acoef[ilev].array(mfi);
            auto bcoeffab = bcoef[ilev].array(mfi);

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_bcoef(i,j,k,bcoeffab,prob_lo,prob_hi,dx);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                actual_init_abeclap_in(i,j,k,rhsfab,exactfab,acoeffab,bcoeffab,
                                       a,b,prob_lo,prob_hi,dx);
            });

            // For inhomogeneous Neumann, we store dphi/d[xyz] in the ghost cells of
            // a cell-centered MultiFab, even though the data represent face values.

            if (bx.smallEnd(0) == domain.smallEnd(0)) {
                Box const& bxlo = amrex::adjCellLo(bx, 0);
                amrex::ParallelFor(bxlo,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dx_lo(i,j,k,solnfab,prob_lo,dx);
                });
            }

            if (bx.bigEnd(0) == domain.bigEnd(0)) {
                Box const& bxhi = amrex::adjCellHi(bx, 0);
                amrex::ParallelFor(bxhi,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dx_hi(i,j,k,solnfab,prob_lo,dx);
                });
            }
            
            if (bx.smallEnd(1) == domain.smallEnd(1)) {
                Box const& bylo = amrex::adjCellLo(bx, 1);
                amrex::ParallelFor(bylo,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dy_lo(i,j,k,solnfab,prob_lo,dx);
                });
            }

            if (bx.bigEnd(1) == domain.bigEnd(1)) {
                Box const& byhi = amrex::adjCellHi(bx, 1);
                amrex::ParallelFor(byhi,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dy_hi(i,j,k,solnfab,prob_lo,dx);
                });
            }

            if (bx.smallEnd(2) == domain.smallEnd(2)) {
                Box const& bzlo = amrex::adjCellLo(bx, 2);
                amrex::ParallelFor(bzlo,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dz_lo(i,j,k,solnfab,prob_lo,dx);
                });
            }

            if (bx.bigEnd(2) == domain.bigEnd(2)) {
                Box const& bzhi = amrex::adjCellHi(bx, 2);
                amrex::ParallelFor(bzhi,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    actual_init_dphi_dz_hi(i,j,k,solnfab,prob_lo,dx);
                });
            }
        }

        solution[ilev].setVal(0.0,0,1,0); // set interior to 0
    }
}
