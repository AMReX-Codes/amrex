#ifndef AMREX_MULTIFAB_UTIL_1D_C_H_
#define AMREX_MULTIFAB_UTIL_1D_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_nd_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& ndfab,
                         int cccomp, int ndcomp, int ncomp) noexcept
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
inline
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& exfab, int cccomp) noexcept
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
inline
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, int cccomp, GeometryData const& gd) noexcept
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

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, FArrayBox& fxfab,
                         FArrayBox const& ccfab, GeometryData const& gd) noexcept
{
    const auto lo = lbound(ndbx);
    const auto fx = fxfab.view(lo);
    const auto cc = ccfab.view(lo);

    const auto xlen = length(ndbx,xbx);

    const int coord_type = gd.Coord();

    switch (coord_type)
    {
    case 0:
    {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < xlen.x; ++i) {
            fx(i,0,0) = 0.5*(cc(i-1,0,0) + cc(i,0,0));
        }
        break;
    }
    case 1:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < xlen.x; ++i) {
            Real rlo = problo + (i+lo.x-0.5)*dx;
            Real rhi = problo + (i+lo.x+0.5)*dx;
            Real rcen = 0.5*(rlo+rhi);
            fx(i,0,0) = 0.5*(rlo*cc(i-1,0,0) + rhi*cc(i,0,0)) / rcen;
        }
        break;
    }
    case 2:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < xlen.x; ++i) {
            Real rlo = problo + (i+lo.x-0.5)*dx;
            Real rhi = problo + (i+lo.x+0.5)*dx;
            Real rcen = 0.5*(rlo+rhi);
            fx(i,0,0) = 0.5*(rlo*rlo*cc(i-1,0,0) + rhi*rhi*cc(i,0,0)) / (rcen*rcen);
        }
        break;
    }
    default:
        amrex::Abort("amrex_avg_cc_to_fc: wrong coord_type");
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_faces (Box const& bx, FArrayBox& crsefab, FArrayBox const& finefab,
                          int ccomp, int fcomp, int ncomp, IntVect const& ratio, int idir) noexcept
{
    const auto len = length(bx);
    const auto clo = lbound(bx);
    const auto flo = refine(clo,ratio);
    const auto crse = crsefab.view(clo,ccomp);
    const auto fine = finefab.view(flo,fcomp);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        for (int i = 0; i < len.x; ++i) {
            crse(i,0,0,n) = fine(facx*i,0,0,n);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_edges (Box const& bx, FArrayBox& crsefab, FArrayBox const& finefab,
                          int ccomp, int fcomp, int ncomp, IntVect const& ratio, int idir) noexcept
{
    const auto len = length(bx);
    const auto clo = lbound(bx);
    const auto flo = refine(clo,ratio);
    const auto crse = crsefab.view(clo,ccomp);
    const auto fine = finefab.view(flo,fcomp);

    const int facx = ratio[0];
    Real facInv = 1./facx;

    for (int n = 0; n < ncomp; ++n) {
        for (int i = 0; i < len.x; ++i) {
            Real c = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(facx*i+iref,0,0,n);
            }
            crse(i,0,0,n) = c * facInv;
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown (Box const& bx, FArrayBox& crsefab, FArrayBox const& finefab,
                    int ccomp, int fcomp, int ncomp, IntVect const& ratio) noexcept
{
    const auto len = length(bx);
    const auto clo = lbound(bx);
    const auto flo = refine(clo,ratio);
    const auto crse = crsefab.view(clo,ccomp);
    const auto fine = finefab.view(flo,fcomp);

    const int facx = ratio[0];
    const Real volfrac = 1.0/static_cast<Real>(facx);

    for (int n = 0; n < ncomp; ++n) {
        for (int i = 0; i < len.x; ++i) {
            int ii = i*facx;
            Real c = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,0,0,n);
            }
            crse(i,0,0,n) = volfrac * c;
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_with_vol (Box const& bx, FArrayBox& crsefab, FArrayBox const& finefab,
                             FArrayBox const& finevolfab, int ccomp, int fcomp, int ncomp,
                             IntVect const& ratio) noexcept
{
    const auto len = length(bx);
    const auto clo = lbound(bx);
    const auto flo = refine(clo,ratio);
    const auto crse = crsefab.view(clo,ccomp);
    const auto fine = finefab.view(flo,fcomp);
    const auto fv = finevolfab.view(flo);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        for (int i = 0; i < len.x; ++i) {
            int ii = i*facx;
            Real cd = 0., cv = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                cv +=                     fv(ii+iref,0,0);
                cd += fine(ii+iref,0,0,n)*fv(ii+iref,0,0);
            }
            crse(i,0,0,n) = cd/cv;
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_nodes (Box const& bx, FArrayBox& crsefab, FArrayBox const& finefab,
                          int ccomp, int fcomp, int ncomp, IntVect const& ratio) noexcept
{
    const auto len = length(bx);
    const auto clo = lbound(bx);
    const auto flo = refine(clo,ratio);
    const auto crse = crsefab.view(clo,ccomp);
    const auto fine = finefab.view(flo,fcomp);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            crse(i,0,0,n) = fine(i*facx,0,0,n);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, FArrayBox& divufab, FArrayBox const& ufab,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto divu = divufab.view(lo);
    const auto    u =    ufab.view(lo);

    const Real dxi = dxinv[0];

    AMREX_PRAGMA_SIMD
    for (int i = 0; i < len.x; ++i) {
        divu(i,0,0) = dxi * (u(i+1,0,0)-u(i,0,0));
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient (Box const& bx, Array4<Real> const& grad,
                             Array4<Real const> const& u,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        grad(i,0,0) = dxi * (u(i+1,0,0)-u(i,0,0));
    }
}

}

#endif
