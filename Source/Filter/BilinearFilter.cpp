#include <WarpX.H>
#include <BilinearFilter.H>
#include <WarpX_f.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

namespace {
    void compute_stencil(Gpu::ManagedVector<Real> &stencil, int npass)
    {
        Gpu::ManagedVector<Real> old_s(1+npass,0.);
        Gpu::ManagedVector<Real> new_s(1+npass,0.);

        old_s[0] = 1.;
        int jmax = 1;
        amrex::Real loc;
        // Convolve the filter with itself npass times
        for(int ipass=1; ipass<npass+1; ipass++){
            // element 0 has to be treated in its own way
            new_s[0] = 0.5 * old_s[0];
            if (1<jmax) new_s[0] += 0.5 * old_s[1];
            loc = 0.;
            // For each element j, apply the filter to 
            // old_s to get new_s[j]. loc stores the tmp 
            // filtered value.
            for(int j=1; j<jmax+1; j++){
                loc = 0.5 * old_s[j];
                loc += 0.25 * old_s[j-1];
                if (j<jmax) loc += 0.25 * old_s[j+1];
                new_s[j] = loc;
            }
            // copy new_s into old_s
            old_s = new_s;
            // extend the stencil length for next iteration
            jmax += 1;
        }
        // we use old_s here to make sure the stencil
        // is corrent even when npass = 0
        stencil = old_s;
        stencil[0] *= 0.5; // because we will use it twice
    }
}

void BilinearFilter::ComputeStencils(){
    BL_PROFILE("BilinearFilter::ComputeStencils()");
    stencil_length_each_dir = npass_each_dir;
    stencil_length_each_dir += 1.;
#if (AMREX_SPACEDIM == 3)
    // npass_each_dir = npass_x npass_y npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_y.resize( 1 + npass_each_dir[1] );
    stencil_z.resize( 1 + npass_each_dir[2] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_y, npass_each_dir[1]);
    compute_stencil(stencil_z, npass_each_dir[2]);
#elif (AMREX_SPACEDIM == 2)
    // npass_each_dir = npass_x npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_z.resize( 1 + npass_each_dir[1] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_z, npass_each_dir[1]);
#endif
    slen = stencil_length_each_dir.dim3();
#if (AMREX_SPACEDIM == 2)
    slen.z = 1;
#endif
}


void
BilinearFilter::ApplyStencil (MultiFab& dstmf, const MultiFab& srcmf, int scomp, int dcomp, int ncomp)
{
    BL_PROFILE("BilinearFilter::ApplyStencil()");
    ncomp = std::min(ncomp, srcmf.nComp());
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox tmpfab;
        for (MFIter mfi(dstmf,TilingIfNotGPU()); mfi.isValid(); ++mfi){
            const auto& srcfab = srcmf[mfi];
            auto& dstfab = dstmf[mfi];
            const Box& tbx = mfi.growntilebox();
            const Box& gbx = amrex::grow(tbx,stencil_length_each_dir-1);
            // tmpfab has enough ghost cells for the stencil
            AsyncFab tmp_async_fab(tmpfab,gbx,ncomp);
            FArrayBox* tmpfab_ptr = tmp_async_fab.fabPtr();
            const FArrayBox* srcfab_ptr = srcmf.fabPtr(mfi);
            // Copy values in srcfab into tmpfab
            const Box& ibx = gbx & srcfab.box();
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(gbx, tgbx,
            {
                tmpfab_ptr->setVal(0.0, tgbx, 0, ncomp);
            });

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(ibx, tibx,
            {
                tmpfab_ptr->copy(*srcfab_ptr, tibx, scomp, tibx, 0, ncomp);
            });

            // Apply filter
            Filter(tbx, tmp_async_fab.hostFab(), dstfab, 0, dcomp, ncomp);
        }
    }
}

void BilinearFilter::Filter (const Box& tbx, FArrayBox const& tmpfab, FArrayBox &dstfab,
                             int scomp, int dcomp, int ncomp)
{
    const auto lo = amrex::lbound(tbx);
    const auto hi = amrex::ubound(tbx);
    const auto tmp = tmpfab.array();
    const auto dst = dstfab.array();
    // tmp and dst are of type Array4 (Fortran ordering)
    amrex::Real const* AMREX_RESTRICT sx = stencil_x.dataPtr();
    amrex::Real const* AMREX_RESTRICT sy = stencil_y.dataPtr();
    amrex::Real const* AMREX_RESTRICT sz = stencil_z.dataPtr();
#ifdef AMREX_USE_CUDA
    Dim3 slen_local = slen;
    amrex::ParallelFor(tbx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dst(i,j,k,dcomp+n) = 0.0;

        for         (int iz=0; iz < slen_local.z; ++iz){
            for     (int iy=0; iy < slen_local.y; ++iy){
                for (int ix=0; ix < slen_local.x; ++ix){
#if (AMREX_SPACEDIM == 3)        
                    Real sss = sx[ix]*sy[iy]*sz[iz];
#else
                    Real sss = sx[ix]*sz[iy];
#endif                        
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
    });
#else  // if not USE_CUDA
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
#endif // USE_CUDA
}
