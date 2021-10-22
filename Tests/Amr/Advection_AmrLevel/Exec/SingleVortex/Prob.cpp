
#include <AMReX_MultiFab.H>

/**
 * Initalize problem
 *
 * \param S_tmp
 * \param geom
 *
 */

using namespace amrex;

//AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void initdata (MultiFab& S_tmp, const Geometry& geom){


    int dm = AMREX_SPACEDIM; //TODO

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();

    for (MFIter mfi(S_tmp); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const Array4<Real>& phi = S_tmp.array(mfi);

        // Populate MultiFab data
        ParallelFor(box, [=] AMREX_GPU_DEVICE ( int i, int j, int k) noexcept
        {
        // AMREX_SPACEDIM
            Real x,y,z,r2;

            x = prob_lo[0] + (i + 0.5) * dx[0];
            y = prob_lo[1] + (j + 0.5) * dx[1];
            z = prob_lo[2] + (k + 0.5) * dx[2];

            if ( dm == 2 ) {
                r2 = ((x-0.5)*(x-0.5) + (y-0.75)*(y-0.75)) / 0.01;
                phi(i,j,k) = 1.0 + std::exp(-r2);
            } else {
               r2 = ((x-0.5)*(x-0.5) + (y-0.75)*(y-0.75) + (z-0.5)*(z-0.5)) / 0.01;
               phi(i,j,k) = 1.0 + std::exp(-r2);
            }

        });
    }

}

