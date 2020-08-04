
#include <AMReX_Extrapolater.H>
#include <AMReX_extrapolater_K.H>
#include <AMReX_iMultiFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace Extrapolater
{
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp)
    {
        BL_ASSERT(mf.nGrow() == 1);
        BL_ASSERT(scomp >= 0);
        BL_ASSERT((scomp+ncomp) <= mf.nComp());

        iMultiFab mask(mf.boxArray(), mf.DistributionMap(), 1, 1, MFInfo(),
                       DefaultFabFactory<IArrayBox>());
        mask.BuildMask(geom.Domain(), geom.periodicity(),
                       finebnd, crsebnd, physbnd, interior);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
           const Box& bx        = mfi.validbox();
           auto const& mask_arr = mask.const_array(mfi);
           auto const& data_arr = mf.array(mfi,scomp);

           if (Gpu::inLaunchRegion()) {
              ParallelFor(amrex::grow(bx,1), ncomp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                 if (mask_arr(i,j,k) == crsebnd) data_arr(i,j,k,n) = 0.0;
              });
              ParallelFor(amrex::grow(bx,1), ncomp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                 amrex_first_order_extrap_gpu(i, j, k, n, bx, mask_arr, data_arr);
              });
           } else {
              amrex_first_order_extrap_cpu(bx, ncomp, mask_arr, data_arr);
           }
        }
    }
}

}
