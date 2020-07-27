
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
            const Box& bx      = mfi.validbox();
            auto const& mask_arr = mask.array(mfi);
            auto const& data_arr = mf.array(mfi,scomp);

            amrex::launch(bx, [mask_arr,data_arr,ncomp]
            AMREX_GPU_DEVICE (Box const& tbx)
            {
               amrex_first_order_extrap(tbx, ncomp, mask_arr, data_arr);
            });
        }
    }
}

}
