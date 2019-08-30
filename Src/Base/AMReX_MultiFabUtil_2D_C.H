#ifndef AMREX_MULTIFAB_UTIL_2D_C_H_
#define AMREX_MULTIFAB_UTIL_2D_C_H_

#include <AMReX_Gpu.H>
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
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            cc(i,j,0,n+cccomp) = 0.25*( nd(i,j  ,0,n+ndcomp) + nd(i+1,j  ,0,n+ndcomp)
                                      + nd(i,j+1,0,n+ndcomp) + nd(i+1,j+1,0,n+ndcomp));
        }}
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_eg_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         Array4<Real const> const& Ey,
                         int cccomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        cc(i,j,0,0+cccomp) = 0.5 * ( Ex(i,j,0) + Ex(i,j+1,0) );
        cc(i,j,0,1+cccomp) = 0.5 * ( Ey(i,j,0) + Ey(i+1,j,0) );
    }}
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_fc_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& fx,
                         Array4<Real const> const& fy,
                         int cccomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        cc(i,j,0,0+cccomp) = 0.5 * ( fx(i,j,0) + fx(i+1,j,0) );
        cc(i,j,0,1+cccomp) = 0.5 * ( fy(i,j,0) + fy(i,j+1,0) );
    }}
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, Box const& ybx,
                         Array4<Real> const& fx, Array4<Real> const& fy,
                         Array4<Real const> const& cc) noexcept
{
    const auto ndlo = lbound(ndbx);
    const auto ndhi = ubound(ndbx);
    const auto xlo = amrex::elemwiseMax(ndlo, lbound(xbx));
    const auto ylo = amrex::elemwiseMax(ndlo, lbound(ybx));
    const auto xhi = amrex::elemwiseMin(ndhi, ubound(xbx));
    const auto yhi = amrex::elemwiseMin(ndhi, ubound(ybx));

    for     (int j = xlo.y; j <= xhi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x; i <= xhi.x; ++i) {
            fx(i,j,0) = 0.5*(cc(i-1,j,0) + cc(i,j,0));
        }
    }

    for     (int j = ylo.y; j <= yhi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x; i <= yhi.x; ++i) {
            fy(i,j,0) = 0.5*(cc(i,j-1,0) + cc(i,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_faces (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int idir) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facy;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int jref = 0; jref < facy; ++jref) {
                    c += fine(ii,jj+jref,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / facx;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_edges (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio, int idir) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facx;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / facy;
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for (int jref = 0; jref < facx; ++jref) {
                    c += fine(ii,jj+jref,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
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
    const int facy = ratio[1];
    const Real volfrac = 1.0/static_cast<Real>(facx*facy);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            Real c = 0.;
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,0,n+fcomp);
            }}
            crse(i,j,0,n+ccomp) = volfrac * c;
        }}
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avgdown_nodes (Box const& bx, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, int ncomp,
                          IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = clo.y; j <= chi.y; ++j) {
            int jj = j*facy;
            AMREX_PRAGMA_SIMD
            for (int i = clo.x; i <= chi.x; ++i) {
                crse(i,j,0,n+ccomp) = fine(i*facx,jj,0,n+fcomp);
            }
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
    const int facy = ratio[1];

    for (int n = 0; n < ncomp; ++n) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            Real cd = 0., cv = 0.;
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                cv +=                                 fv(ii+iref,jj+jref,0);
                cd += fine(ii+iref,jj+jref,0,n+fcomp)*fv(ii+iref,jj+jref,0);
            }}
            crse(i,j,0,n+ccomp) = cd/cv;
        }}
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, Array4<Real> const& divu,
                               Array4<Real const> const& u,
                               Array4<Real const> const& v,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            divu(i,j,0) = dxi * (u(i+1,j,0)-u(i,j,0))
                +         dyi * (v(i,j+1,0)-v(i,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient (Box const& bx, Array4<Real> const& grad,
                             Array4<Real const> const& u,
                             Array4<Real const> const& v,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            grad(i,j,0,0) = dxi * (u(i+1,j,0)-u(i,j,0));
            grad(i,j,0,1) = dyi * (v(i,j+1,0)-v(i,j,0));
        }
    }
}

}

#endif
