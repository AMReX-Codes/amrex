
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
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,0,0,n) = 0.5*(nd(i,0,0,n)+nd(i+1,0,0,n));
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& exfab, int cccomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto Ex = exfab.view(lo);
    
    AMREX_PRAGMA_SIMD
    for (int i = 0; i < len.x; ++i) {
        cc(i,0,0) = Ex(i,0,0);
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, int cccomp, GeometryData const& gd)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto fx = fxfab.view(lo);

    const int coord_type = gd.Coord();

    switch (coord_type)
    {
    case 0:
    {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,0,0) = 0.5 * ( fx(i,0,0) + fx(i+1,0,0) );
        }
        break;
    }
    case 1:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real rlo = problo + (i+lo.x)*dx;
            Real rhi = problo + (i+1+lo.x)*dx;
            Real rcen = 0.5*(rlo+rhi);
            cc(i,0,0) = 0.5 * ( rlo*fx(i,0,0) + rhi*fx(i+1,0,0) ) / rcen;
        }
        break;
    }
    case 2:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real rlo = problo + (i+lo.x)*dx;
            Real rhi = problo + (i+1+lo.x)*dx;
            Real rcen = 0.5*(rlo+rhi);
            cc(i,0,0) = 0.5 * ( rlo*rlo*fx(i,0,0) + rhi*rhi*fx(i+1,0,0) ) / (rcen*rcen);
        }
        break;
    }
    default:
        amrex::Abort("amrex_avg_fc_to_cc: wrong coord_type");
    }
}

}

