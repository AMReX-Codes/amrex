#include <WarpX.H>
#include <Filter.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

#ifdef AMREX_USE_CUDA

/* \brief Apply stencil on MultiFab (GPU version, 2D/3D).
 * \param dstmf Destination MultiFab
 * \param srcmf source MultiFab
 * \param scomp first component of srcmf on which the filter is applied
 * \param dcomp first component of dstmf on which the filter is applied
 * \param ncomp Number of components on which the filter is applied.
 */
void
Filter::ApplyStencil (MultiFab& dstmf, const MultiFab& srcmf, int scomp, int dcomp, int ncomp)
{
    BL_PROFILE("BilinearFilter::ApplyStencil(MultiFab)");
    ncomp = std::min(ncomp, srcmf.nComp());

    for (MFIter mfi(dstmf); mfi.isValid(); ++mfi)
    {
        const auto& src = srcmf.array(mfi);
        const auto& dst = dstmf.array(mfi);
        const Box& tbx = mfi.growntilebox();
        const Box& gbx = amrex::grow(tbx,stencil_length_each_dir-1);

        // tmpfab has enough ghost cells for the stencil
        FArrayBox tmp_fab(gbx,ncomp);
        Elixir tmp_eli = tmp_fab.elixir();  // Prevent the tmp data from being deleted too early
        auto const& tmp = tmp_fab.array();

        // Copy values in srcfab into tmpfab
        const Box& ibx = gbx & srcmf[mfi].box();
        AMREX_PARALLEL_FOR_4D ( gbx, ncomp, i, j, k, n,
        {
            if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                tmp(i,j,k,n) = src(i,j,k,n+scomp);
            } else {
                tmp(i,j,k,n) = 0.0;
            }
        });

        // Apply filter
        DoFilter(tbx, tmp, dst, 0, dcomp, ncomp);
    }
}

/* \brief Apply stencil on FArrayBox (GPU version, 2D/3D).
 * \param dstfab Destination FArrayBox
 * \param srcmf source FArrayBox
 * \param tbx Grown box on which srcfab is defined.
 * \param scomp first component of srcfab on which the filter is applied
 * \param dcomp first component of dstfab on which the filter is applied
 * \param ncomp Number of components on which the filter is applied.
 */
void
Filter::ApplyStencil (FArrayBox& dstfab, const FArrayBox& srcfab, 
                      const Box& tbx, int scomp, int dcomp, int ncomp)
{
    BL_PROFILE("BilinearFilter::ApplyStencil(FArrayBox)");
    ncomp = std::min(ncomp, srcfab.nComp());
    const auto& src = srcfab.array();
    const auto& dst = dstfab.array();
    const Box& gbx = amrex::grow(tbx,stencil_length_each_dir-1);

    // tmpfab has enough ghost cells for the stencil
    FArrayBox tmp_fab(gbx,ncomp);
    Elixir tmp_eli = tmp_fab.elixir();  // Prevent the tmp data from being deleted too early
    auto const& tmp = tmp_fab.array();

    // Copy values in srcfab into tmpfab
    const Box& ibx = gbx & srcfab.box();
    AMREX_PARALLEL_FOR_4D ( gbx, ncomp, i, j, k, n,
        {
            if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                tmp(i,j,k,n) = src(i,j,k,n+scomp);
            } else {
                tmp(i,j,k,n) = 0.0;
            }
        });

    // Apply filter
    DoFilter(tbx, tmp, dst, 0, dcomp, ncomp);
}

/* \brief Apply stencil (2D/3D, CPU/GPU)
 */
void Filter::DoFilter (const Box& tbx,
                       Array4<Real const> const& tmp,
                       Array4<Real      > const& dst,
                       int scomp, int dcomp, int ncomp)
{
    amrex::Real const* AMREX_RESTRICT sx = stencil_x.data();
    amrex::Real const* AMREX_RESTRICT sy = stencil_y.data();
    amrex::Real const* AMREX_RESTRICT sz = stencil_z.data();
    Dim3 slen_local = slen;
    AMREX_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
    {
        Real d = 0.0;

        for         (int iz=0; iz < slen_local.z; ++iz){
            for     (int iy=0; iy < slen_local.y; ++iy){
                for (int ix=0; ix < slen_local.x; ++ix){
#if (AMREX_SPACEDIM == 3)        
                    Real sss = sx[ix]*sy[iy]*sz[iz];
#else
                    Real sss = sx[ix]*sz[iy];
#endif                        
#if (AMREX_SPACEDIM == 3)
                    d += sss*( tmp(i-ix,j-iy,k-iz,scomp+n)
                              +tmp(i+ix,j-iy,k-iz,scomp+n)
                              +tmp(i-ix,j+iy,k-iz,scomp+n)
                              +tmp(i+ix,j+iy,k-iz,scomp+n)
                              +tmp(i-ix,j-iy,k+iz,scomp+n)
                              +tmp(i+ix,j-iy,k+iz,scomp+n)
                              +tmp(i-ix,j+iy,k+iz,scomp+n)
                              +tmp(i+ix,j+iy,k+iz,scomp+n));
#else
                    d += sss*( tmp(i-ix,j-iy,k,scomp+n)
                              +tmp(i+ix,j-iy,k,scomp+n)
                              +tmp(i-ix,j+iy,k,scomp+n)
                              +tmp(i+ix,j+iy,k,scomp+n));
#endif
                }
            }
        }

        dst(i,j,k,dcomp+n) = d;
    });
}

#else

/* \brief Apply stencil on MultiFab (CPU version, 2D/3D).
 * \param dstmf Destination MultiFab
 * \param srcmf source MultiFab
 * \param scomp first component of srcmf on which the filter is applied
 * \param dcomp first component of dstmf on which the filter is applied
 * \param ncomp Number of components on which the filter is applied.
 */
void
Filter::ApplyStencil (MultiFab& dstmf, const MultiFab& srcmf, int scomp, int dcomp, int ncomp)
{
    BL_PROFILE("BilinearFilter::ApplyStencil()");
    ncomp = std::min(ncomp, srcmf.nComp());
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox tmpfab;
        for (MFIter mfi(dstmf,true); mfi.isValid(); ++mfi){
            const auto& srcfab = srcmf[mfi];
            auto& dstfab = dstmf[mfi];
            const Box& tbx = mfi.growntilebox();
            const Box& gbx = amrex::grow(tbx,stencil_length_each_dir-1);
            // tmpfab has enough ghost cells for the stencil
            tmpfab.resize(gbx,ncomp);
            tmpfab.setVal(0.0, gbx, 0, ncomp);
            // Copy values in srcfab into tmpfab
            const Box& ibx = gbx & srcfab.box();
            tmpfab.copy(srcfab, ibx, scomp, ibx, 0, ncomp);
            // Apply filter
            DoFilter(tbx, tmpfab.array(), dstfab.array(), 0, dcomp, ncomp);
        }
    }
}

/* \brief Apply stencil on FArrayBox (CPU version, 2D/3D).
 * \param dstfab Destination FArrayBox
 * \param srcmf source FArrayBox
 * \param tbx Grown box on which srcfab is defined.
 * \param scomp first component of srcfab on which the filter is applied
 * \param dcomp first component of dstfab on which the filter is applied
 * \param ncomp Number of components on which the filter is applied.
 */
void
Filter::ApplyStencil (FArrayBox& dstfab, const FArrayBox& srcfab, 
                      const Box& tbx, int scomp, int dcomp, int ncomp)
{
    BL_PROFILE("BilinearFilter::ApplyStencil(FArrayBox)");
    ncomp = std::min(ncomp, srcfab.nComp());
    FArrayBox tmpfab;
    const Box& gbx = amrex::grow(tbx,stencil_length_each_dir-1);
    // tmpfab has enough ghost cells for the stencil
    tmpfab.resize(gbx,ncomp);
    tmpfab.setVal(0.0, gbx, 0, ncomp);
    // Copy values in srcfab into tmpfab
    const Box& ibx = gbx & srcfab.box();
    tmpfab.copy(srcfab, ibx, scomp, ibx, 0, ncomp);
    // Apply filter
    DoFilter(tbx, tmpfab.array(), dstfab.array(), 0, dcomp, ncomp);
}

void Filter::DoFilter (const Box& tbx,
                       Array4<Real const> const& tmp,
                       Array4<Real      > const& dst,
                       int scomp, int dcomp, int ncomp)
{
    const auto lo = amrex::lbound(tbx);
    const auto hi = amrex::ubound(tbx);
    // tmp and dst are of type Array4 (Fortran ordering)
    amrex::Real const* AMREX_RESTRICT sx = stencil_x.data();
    amrex::Real const* AMREX_RESTRICT sy = stencil_y.data();
    amrex::Real const* AMREX_RESTRICT sz = stencil_z.data();
    for (int n = 0; n < ncomp; ++n) {
        // Set dst value to 0.
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    dst(i,j,k,dcomp+n) = 0.0;
                }
            }
        }
        // 3 nested loop on 3D stencil
        for         (int iz=0; iz < slen.z; ++iz){
            for     (int iy=0; iy < slen.y; ++iy){
                for (int ix=0; ix < slen.x; ++ix){
#if (AMREX_SPACEDIM == 3)        
                    Real sss = sx[ix]*sy[iy]*sz[iz];
#else
                    Real sss = sx[ix]*sz[iy];
#endif
                    // 3 nested loop on 3D array
                    for         (int k = lo.z; k <= hi.z; ++k) {
                        for     (int j = lo.y; j <= hi.y; ++j) {
                            AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {
#if (AMREX_SPACEDIM == 3)
                                dst(i,j,k,dcomp+n) += sss*(tmp(i-ix,j-iy,k-iz,scomp+n)
                                                          +tmp(i+ix,j-iy,k-iz,scomp+n)
                                                          +tmp(i-ix,j+iy,k-iz,scomp+n)
                                                          +tmp(i+ix,j+iy,k-iz,scomp+n)
                                                          +tmp(i-ix,j-iy,k+iz,scomp+n)
                                                          +tmp(i+ix,j-iy,k+iz,scomp+n)
                                                          +tmp(i-ix,j+iy,k+iz,scomp+n)
                                                          +tmp(i+ix,j+iy,k+iz,scomp+n));
#else
                                dst(i,j,k,dcomp+n) += sss*(tmp(i-ix,j-iy,k,scomp+n)
                                                          +tmp(i+ix,j-iy,k,scomp+n)
                                                          +tmp(i-ix,j+iy,k,scomp+n)
                                                          +tmp(i+ix,j+iy,k,scomp+n));
#endif
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif // #ifdef AMREX_USE_CUDA
