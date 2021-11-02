
#include <AMReX_MultiFab.H>

/**
 * Initialize Data on Multifab
 *
 * \param S_tmp Pointer to Multifab where data is to be initialized.
 * \param geom Pointer to Multifab's geometry data.
 *
 */

using namespace amrex;

void initdata (MultiFab& S_tmp, const Geometry& geom){

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();

    for (MFIter mfi(S_tmp); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const Array4<Real>& phi = S_tmp.array(mfi);

        // Populate MultiFab data
        ParallelFor(box, [=] AMREX_GPU_DEVICE ( int i, int j, int k) noexcept
        {

            Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];

#if (AMREX_SPACEDIM == 2)

            Real r2 = ((x-0.5_rt)*(x-0.5_rt) + (y-0.75_rt)*(y-0.75_rt)) / 0.01_rt;
            phi(i,j,k) = 1.0_rt + std::exp(-r2);

#elif (AMREX_SPACEDIM == 3)

            Real z = prob_lo[2] + (k + 0.5_rt) * dx[2];
            Real r2 = ((x-0.5_rt)*(x-0.5_rt) + (y-0.75_rt)*(y-0.75_rt) + (z-0.5_rt)*(z-0.5_rt)) / 0.01_rt;
            phi(i,j,k) = 1.0_rt + std::exp(-r2);
#endif

        });
    }

}

