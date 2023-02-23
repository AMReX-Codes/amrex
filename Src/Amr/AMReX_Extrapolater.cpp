
#include <AMReX_Extrapolater.H>
#include <AMReX_extrapolater_K.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex::Extrapolater
{
    // Backward compatible version filling only 1 ghost cell
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp)
    {
       FirstOrderExtrap(mf, geom, scomp, ncomp, 1);
    }

    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp, int ngrow)
    {
        BL_ASSERT(mf.nGrow() >= ngrow);
        BL_ASSERT(scomp >= 0);
        BL_ASSERT((scomp+ncomp) <= mf.nComp());

        iMultiFab mask(mf.boxArray(), mf.DistributionMap(), 1, ngrow, MFInfo(),
                       DefaultFabFactory<IArrayBox>());
        mask.BuildMask(geom.Domain(), geom.periodicity(),
                       finebnd, crsebnd, physbnd, interior);

        // Do the extrap. on successive layers of ghost cells
        for (int layer = 0; layer < ngrow; layer++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Box& bx        = mfi.validbox();
                const Box& gbx       = amrex::grow(bx,layer);
                auto const& mask_arr = mask.const_array(mfi);
                auto const& data_arr = mf.array(mfi,scomp);

                if (Gpu::inLaunchRegion()) {
                    // set the crse cell to zero in the current layer
                    ParallelFor(amrex::grow(gbx,1), ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                       if (!gbx.contains(i,j,k)) {
                           if (mask_arr(i,j,k) == crsebnd) data_arr(i,j,k,n) = 0.0;
                       }
                    });
                    ParallelFor(amrex::grow(gbx,1), ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                       amrex_first_order_extrap_gpu(i, j, k, n, bx, mask_arr, data_arr);
                    });
                } else {
                    amrex_first_order_extrap_cpu(gbx, ncomp, mask_arr, data_arr);
                }
            }
        }
    }
}
