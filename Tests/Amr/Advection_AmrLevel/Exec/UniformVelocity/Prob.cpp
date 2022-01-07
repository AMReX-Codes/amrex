
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

            Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
            Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

#if (AMREX_SPACEDIM == 2)

            Real r2 = ((x-Real(0.0))*(x-Real(0.0)) + (y-Real(0.0))*(y-Real(0.0))) / Real(0.01);
            phi(i,j,k) = Real(1.0) + std::exp(-r2);

#elif (AMREX_SPACEDIM == 3)

            Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
            Real r2 = ((x-Real(0.0))*(x-Real(0.0)) + (y-Real(0.0))*(y-Real(0.0)) + (z-Real(0.0))*(z-Real(0.0))) / Real(0.01);
            phi(i,j,k) = Real(1.0) + std::exp(-r2);
#endif

        });
    }

}

