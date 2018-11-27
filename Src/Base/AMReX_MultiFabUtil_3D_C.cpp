
#include <AMReX_MultiFabUtil_C.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_avg_nd_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& ndfab,
                         int cccomp, int ndcomp, int ncomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto nd = ndfab.view(lo,ndcomp);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
        for (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,k,n) = 0.125*( nd(i,j  ,k  ,n) + nd(i+1,j  ,k  ,n)
                                + nd(i,j+1,k  ,n) + nd(i+1,j+1,k  ,n)
                                + nd(i,j  ,k+1,n) + nd(i+1,j  ,k+1,n)
                                + nd(i,j+1,k+1,n) + nd(i+1,j+1,k+1,n));
        }}}
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& exfab, FArrayBox const& eyfab, FArrayBox const& ezfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo);
    const auto Ex = exfab.view(lo);
    const auto Ey = eyfab.view(lo);
    const auto Ez = ezfab.view(lo);
    
    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                cc(i,j,k,0) = 0.25 * ( Ex(i,j,k) + Ex(i,j+1,k) + Ex(i,j,k+1) + Ex(i,j+1,k+1) );
                cc(i,j,k,1) = 0.25 * ( Ey(i,j,k) + Ey(i+1,j,k) + Ey(i,j,k+1) + Ey(i+1,j,k+1) );
                cc(i,j,k,2) = 0.25 * ( Ez(i,j,k) + Ez(i+1,j,k) + Ez(i,j+1,k) + Ez(i+1,j+1,k) );
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, FArrayBox const& fyfab, FArrayBox const& fzfab,
                         GeometryData const& gd)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo);
    const auto fx = fxfab.view(lo);
    const auto fy = fyfab.view(lo);
    const auto fz = fzfab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                cc(i,j,k,0) = 0.5 * ( fx(i,j,k) + fx(i+1,j,k) );
                cc(i,j,k,1) = 0.5 * ( fy(i,j,k) + fy(i,j+1,k) );
                cc(i,j,k,2) = 0.5 * ( fz(i,j,k) + fz(i,j,k+1) );
            }
        }
    }
}

}
