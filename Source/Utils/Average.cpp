#include "Average.H"

using namespace amrex;

void
Average::CoarsenAndInterpolateLoop ( MultiFab& mf_dst,
                                     const MultiFab& mf_src,
                                     const int dcomp,
                                     const int scomp,
                                     const int ncomp,
                                     const int ngrow,
                                     const IntVect crse_ratio )
{
    // Staggering of source fine MultiFab and destination coarse MultiFab
    const IntVect stag_src = mf_src.boxArray().ixType().toIntVect();
    const IntVect stag_dst = mf_dst.boxArray().ixType().toIntVect();

    if ( crse_ratio > IntVect(1) ) AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ngrow == 0,
        "option of filling guard cells of destination MultiFab with coarsening not supported for this interpolation" );

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( mf_src.nGrowVect() >= stag_dst-stag_src+IntVect(ngrow),
        "source fine MultiFab does not have enough guard cells for this interpolation" );

    // Auxiliary integer arrays (always 3D)
    Gpu::ManagedVector<int> sf_gpuarr, sc_gpuarr, cr_gpuarr, np_gpuarr;
    sf_gpuarr.resize( 3 ); // staggering of source fine MultiFab
    sc_gpuarr.resize( 3 ); // staggering of destination coarse MultiFab
    cr_gpuarr.resize( 3 ); // coarsening ratio
    np_gpuarr.resize( 3 ); // number of points to loop over for interpolation

    sf_gpuarr[0] = stag_src[0];
    sf_gpuarr[1] = stag_src[1];
#if   (AMREX_SPACEDIM == 2)
    sf_gpuarr[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sf_gpuarr[2] = stag_src[2];
#endif

    sc_gpuarr[0] = stag_dst[0];
    sc_gpuarr[1] = stag_dst[1];
#if   (AMREX_SPACEDIM == 2)
    sc_gpuarr[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sc_gpuarr[2] = stag_dst[2];
#endif

    cr_gpuarr[0] = crse_ratio[0];
    cr_gpuarr[1] = crse_ratio[1];
#if   (AMREX_SPACEDIM == 2)
    cr_gpuarr[2] = 1;
#elif (AMREX_SPACEDIM == 3)
    cr_gpuarr[2] = crse_ratio[2];
#endif

    int const* const AMREX_RESTRICT sf = sf_gpuarr.data();
    int const* const AMREX_RESTRICT sc = sc_gpuarr.data();
    int const* const AMREX_RESTRICT cr = cr_gpuarr.data();
    int      * const AMREX_RESTRICT np = np_gpuarr.data();

    // Compute number of points to loop over (either 1 or 2) in each direction
    for ( int l = 0; l < 3; ++l ) {
        if ( cr[l] == 1 ) np[l] = 1+abs(sf[l]-sc[l]); // no coarsening
        else              np[l] = 2-sf[l];
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_dst, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        // Tiles defined at the coarse level
        const Box& bx = mfi.growntilebox( ngrow );
        Array4<Real> const& arr_dst = mf_dst.array( mfi );
        Array4<Real const> const& arr_src = mf_src.const_array( mfi );
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         arr_dst(i,j,k,n+dcomp) = Average::CoarsenAndInterpolateKernel(
                             arr_src, sf, sc, cr, np, i, j, k, n+scomp );
                     } );
    }
}

void
Average::CoarsenAndInterpolate ( MultiFab& mf_dst,
                                 const MultiFab& mf_src,
                                 const int dcomp,
                                 const int scomp,
                                 const int ncomp,
                                 const int ngrow,
                                 const IntVect crse_ratio )
{
    BL_PROFILE( "Average::CoarsenAndInterpolate" );

    // Convert BoxArray of source MultiFab to staggering of destination MultiFab and coarsen it (if possible)
    BoxArray ba_tmp = amrex::convert( mf_src.boxArray(), mf_dst.ixType().toIntVect() );
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ba_tmp.coarsenable( crse_ratio ),
        "source MultiFab converted to staggering of destination MultiFab is not coarsenable" );
    ba_tmp.coarsen( crse_ratio );

    if ( ba_tmp == mf_dst.boxArray() and mf_src.DistributionMap() == mf_dst.DistributionMap() )
        Average::CoarsenAndInterpolateLoop( mf_dst, mf_src, dcomp, scomp, ncomp, ngrow, crse_ratio );
    else
    {
        // Cannot coarsen into MultiFab with different BoxArray or DistributionMapping:
        // 1) create temporary MultiFab on coarsened version of source BoxArray with same DistributionMapping
        MultiFab mf_tmp( ba_tmp, mf_src.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory() );
        // 2) interpolate from mf_src to mf_tmp (start writing into component 0)
        Average::CoarsenAndInterpolateLoop( mf_tmp, mf_src, 0, scomp, ncomp, ngrow, crse_ratio );
        // 3) copy from mf_tmp to mf_dst (with different BoxArray or DistributionMapping)
        mf_dst.copy( mf_tmp, 0, dcomp, ncomp );
    }
}
