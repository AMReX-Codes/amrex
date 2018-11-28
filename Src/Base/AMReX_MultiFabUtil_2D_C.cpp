
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
        for (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,0,n) = 0.25*(nd(i,j,0,n)+nd(i+1,j,0,n)+nd(i,j+1,0,n)+nd(i+1,j+1,0,n));
        }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& exfab, FArrayBox const& eyfab,
                         int cccomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto Ex = exfab.view(lo);
    const auto Ey = eyfab.view(lo);
    
    for     (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,0,0) = 0.5 * ( Ex(i,j,0) + Ex(i,j+1,0) );
            cc(i,j,0,1) = 0.5 * ( Ey(i,j,0) + Ey(i+1,j,0) );
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, FArrayBox const& fyfab,
                         int cccomp, GeometryData const& gd)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto fx = fxfab.view(lo);
    const auto fy = fyfab.view(lo);

    for     (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,0,0) = 0.5 * ( fx(i,j,0) + fx(i+1,j,0) );
            cc(i,j,0,1) = 0.5 * ( fy(i,j,0) + fy(i,j+1,0) );
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, Box const& ybx,
                         FArrayBox& fxfab, FArrayBox& fyfab,
                         FArrayBox const& ccfab, GeometryData const& gd)
{
    const auto lo = lbound(ndbx);
    const auto fx = fxfab.view(lo);
    const auto fy = fyfab.view(lo);
    const auto cc = ccfab.view(lo);

    const auto xlen = length(ndbx,xbx);
    for     (int j = 0; j < xlen.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < xlen.x; ++i) {
            fx(i,j,0) = 0.5*(cc(i-1,j,0) + cc(i,j,0));
        }
    }

    const auto ylen = length(ndbx,ybx);
    for     (int j = 0; j < ylen.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < ylen.x; ++i) {
            fy(i,j,0) = 0.5*(cc(i,j-1,0) + cc(i,j,0));
        }
    }
}

}
