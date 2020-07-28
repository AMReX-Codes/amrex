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
void amrex_avg_nd_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& nd,
                         int cccomp, int ndcomp, int ncomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            cc(i,0,0,n+cccomp) = 0.5*(nd(i,0,0,n+ndcomp)+nd(i+1,0,0,n+ndcomp));
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_eg_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         int cccomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        cc(i,0,0,cccomp) = Ex(i,0,0);
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_fc_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& fx,
                         int cccomp, GeometryData const& gd) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const int coord_type = gd.Coord();

    switch (coord_type)
    {
    case 0:
    {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            cc(i,0,0,cccomp) = 0.5 * ( fx(i,0,0) + fx(i+1,0,0) );
        }
        break;
    }
    case 1:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real rlo = problo + i*dx;
            Real rhi = problo + (i+1)*dx;
            Real rcen = 0.5*(rlo+rhi);
            cc(i,0,0,cccomp) = 0.5 * ( rlo*fx(i,0,0) + rhi*fx(i+1,0,0) ) / rcen;
        }
        break;
    }
    case 2:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real rlo = problo + i*dx;
            Real rhi = problo + (i+1)*dx;
            Real rcen = 0.5*(rlo+rhi);
            cc(i,0,0,cccomp) = 0.5 * ( rlo*rlo*fx(i,0,0) + rhi*rhi*fx(i+1,0,0) ) / (rcen*rcen);
        }
        break;
    }
    default:
        amrex::Abort("amrex_avg_fc_to_cc: wrong coord_type");
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx,
                         Array4<Real> const& fx,
                         Array4<Real const> const& cc,
                         GeometryData const& gd) noexcept
{
    const auto ndlo = lbound(ndbx);
    const auto ndhi = ubound(ndbx);

    const auto xlo = amrex::elemwiseMax(ndlo, lbound(xbx));
    const auto xhi = amrex::elemwiseMin(ndhi, ubound(xbx));

    const int coord_type = gd.Coord();

    switch (coord_type)
    {
    case 0:
    {
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x; i <= xhi.x; ++i) {
            fx(i,0,0) = 0.5*(cc(i-1,0,0) + cc(i,0,0));
        }
        break;
    }
    case 1:
    {
        const Real problo = gd.ProbLo(0);
        const Real dx = gd.CellSize(0);
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x; i <= xhi.x; ++i) {
            Real rlo = problo + (i-0.5)*dx;
            Real rhi = problo + (i+0.5)*dx;
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
        for (int i = xlo.x; i <= xhi.x; ++i) {
            Real rlo = problo + (i-0.5)*dx;
            Real rhi = problo + (i+0.5)*dx;
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
void amrex_avgdown_faces (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int /*idir*/) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        for (int i = clo.x; i <= chi.x; ++i) {
            crse(i,0,0,n+ccomp) = fine(facx*i,0,0,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_edges (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int /*idir*/) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);

    const int facx = ratio[0];
    Real facInv = 1./facx;

    for (int n = 0; n < ncomp; ++n) {
        for (int i = clo.x; i <= chi.x; ++i) {
            Real c = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(facx*i+iref,0,0,n+fcomp);
            }
            crse(i,0,0,n+ccomp) = c * facInv;
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown (Box const& bx, Array4<Real> const& crse,
                    Array4<Real const> const& fine,
                    int ccomp, int fcomp, int ncomp,
                    IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);

    const int facx = ratio[0];
    const Real volfrac = 1.0/static_cast<Real>(facx);

    for (int n = 0; n < ncomp; ++n) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            Real c = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,0,0,n+fcomp);
            }
            crse(i,0,0,n+ccomp) = volfrac * c;
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_with_vol (Box const& bx, Array4<Real> const& crse,
                             Array4<Real const> const& fine,
                             Array4<Real const> const& fv,
                             int ccomp, int fcomp, int ncomp,
                             IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            Real cd = 0., cv = 0.;
            for (int iref = 0; iref < facx; ++iref) {
                cv +=                           fv(ii+iref,0,0);
                cd += fine(ii+iref,0,0,fcomp+n)*fv(ii+iref,0,0);
            }
            crse(i,0,0,ccomp+n) = cd/cv;
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_nodes (Box const& bx, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);

    const int facx = ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = clo.x; i <= chi.x; ++i) {
            crse(i,0,0,n+ccomp) = fine(i*facx,0,0,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, Array4<Real> const& divu,
                               Array4<Real const> const& u,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];

    for     (int n = 0; n < divu.ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            divu(i,0,0,n) = dxi * (u(i+1,0,0,n)-u(i,0,0,n));
        }
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
