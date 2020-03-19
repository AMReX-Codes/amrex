#include "Average.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real Average::ToCellCenter ( Array4<Real const> const& mf_in_arr,
                             const IntVect stag,
                             const int i,
                             const int j,
                             const int k,
                             const int comp)
{
    const int sx = stag[0];
    const int sy = stag[1];
#if   (AMREX_SPACEDIM == 2)
    constexpr int sz = 0;
#elif (AMREX_SPACEDIM == 3)
    const int sz = stag[2];
#endif
    return 0.125_rt * (   mf_in_arr(i   ,j   ,k   ,comp)
                        + mf_in_arr(i+sx,j   ,k   ,comp)
                        + mf_in_arr(i   ,j+sy,k   ,comp)
                        + mf_in_arr(i   ,j   ,k+sz,comp)
                        + mf_in_arr(i+sx,j+sy,k   ,comp)
                        + mf_in_arr(i   ,j+sy,k+sz,comp)
                        + mf_in_arr(i+sx,j   ,k+sz,comp)
                        + mf_in_arr(i+sx,j+sy,k+sz,comp) );
}

void
Average::ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp,
                        const int ncomp )
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_out, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        const Box bx = mfi.growntilebox( ngrow );
        Array4<Real> const& mf_out_arr = mf_out.array( mfi );
        Array4<Real const> const& mf_in_arr = mf_in.const_array( mfi );
        const IntVect stag = mf_in.boxArray().ixType().ixType();
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         mf_out_arr(i,j,k,n+dcomp) = Average::ToCellCenter( mf_in_arr, stag, i, j, k, n+scomp );
                     } );
    }
}
