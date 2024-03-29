#ifndef AMREX_MULTIFAB_UTIL_2D_C_H_
#define AMREX_MULTIFAB_UTIL_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_nd_to_cc (int i, int j, int, int n,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& nd,
                         int cccomp, int ndcomp) noexcept
{
    cc(i,j,0,n+cccomp) = Real(0.25)*( nd(i,j  ,0,n+ndcomp) + nd(i+1,j  ,0,n+ndcomp)
                                    + nd(i,j+1,0,n+ndcomp) + nd(i+1,j+1,0,n+ndcomp));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_eg_to_cc (int i, int j, int,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         Array4<Real const> const& Ey,
                         int cccomp) noexcept
{
    cc(i,j,0,0+cccomp) = Real(0.5) * ( Ex(i,j,0) + Ex(i,j+1,0) );
    cc(i,j,0,1+cccomp) = Real(0.5) * ( Ey(i,j,0) + Ey(i+1,j,0) );
}

template <typename CT, typename FT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_fc_to_cc (int i, int j, int,
                         Array4<CT      > const& cc,
                         Array4<FT const> const& fx,
                         Array4<FT const> const& fy,
                         int cccomp) noexcept
{
    cc(i,j,0,0+cccomp) = CT(0.5) * CT( fx(i,j,0) + fx(i+1,j,0) );
    cc(i,j,0,1+cccomp) = CT(0.5) * CT( fy(i,j,0) + fy(i,j+1,0) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_cc_to_fc (int i, int j, int, int n, Box const& xbx, Box const& ybx,
                         Array4<Real> const& fx, Array4<Real> const& fy,
                         Array4<Real const> const& cc, bool use_harmonic_averaging)  noexcept
{
    if (use_harmonic_averaging)
    {
        if (xbx.contains(i,j,0)) {
            if (cc(i-1,j,0,n) == Real(0.0) || cc(i,j,0,n) == Real(0.0)) {
                fx(i,j,0,n) = Real(0.0);
            } else {
                fx(i,j,0,n) = Real(2.0) / (Real(1.0) / cc(i-1,j,0,n) + Real(1.0) / cc(i,j,0,n));
            }
        }
        if (ybx.contains(i,j,0)) {
            if (cc(i,j-1,0,n) == Real(0.0) || cc(i,j,0,n) == Real(0.0)) {
                fy(i,j,0,n) = Real(0.0);
            } else {
                fy(i,j,0,n) = Real(2.0) / (Real(1.0) / cc(i,j-1,0,n) + Real(1.0) / cc(i,j,0,n));
            }
        }
    } else {
        if (xbx.contains(i,j,0)) {
            fx(i,j,0,n) = Real(0.5)*(cc(i-1,j,0,n) + cc(i,j,0,n));
        }
        if (ybx.contains(i,j,0)) {
            fy(i,j,0,n) = Real(0.5)*(cc(i,j-1,0,n) + cc(i,j,0,n));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_faces (Box const& bx, Array4<T> const& crse,
                          Array4<T const> const& fine,
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
        T facInv =  T(1.0) / static_cast<T>(facy);
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                T c = T(0.);
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
        T facInv = T(1.0) / static_cast<T>(facx);
        for (int n = 0; n < ncomp; ++n) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                T c = T(0.);
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,0,n+fcomp);
                }
                crse(i,j,0,n+ccomp) = c * facInv;
            }}
        }
        break;
    }
    default: { break; }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_faces (int i, int j, int, int n, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        const T facInv =  T(1.0) / static_cast<T>(facy);
        const int ii = i*facx;
        const int jj = j*facy;
        T c = T(0.);
        for (int jref = 0; jref < facy; ++jref) {
            c += fine(ii,jj+jref,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const T facInv = T(1.0) / static_cast<T>(facx);
        const int ii = i*facx;
        const int jj = j*facy;
        T c = T(0.);
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    default: { break; }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
        Real facInv = Real(1.0) / static_cast<Real>(facx);
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
        Real facInv = Real(1.0) / static_cast<Real>(facy);
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
    default: { break; }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_edges (int i, int j, int, int n, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];

    switch (idir) {
    case 0:
    {
        const Real facInv = Real(1.0) / static_cast<Real>(facx);
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const Real facInv = Real(1.0) / static_cast<Real>(facy);
        const int ii = i*facx;
        const int jj = j*facy;
        Real c = 0.;
        for (int jref = 0; jref < facx; ++jref) {
            c += fine(ii,jj+jref,0,n+fcomp);
        }
        crse(i,j,0,n+ccomp) = c * facInv;
        break;
    }
    default: { break; }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown (Box const& bx, Array4<T> const& crse,
                    Array4<T const> const& fine,
                    int ccomp, int fcomp, int ncomp,
                    IntVect const& ratio) noexcept
{
    const auto clo = lbound(bx);
    const auto chi = ubound(bx);
    const int facx = ratio[0];
    const int facy = ratio[1];
    const T volfrac = T(1.0)/T(facx*facy);

    for (int n = 0; n < ncomp; ++n) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            T c = 0;
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,0,n+fcomp);
            }}
            crse(i,j,0,n+ccomp) = volfrac * c;
        }}
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown (int i, int j, int, int n, Array4<T> const& crse,
                    Array4<T const> const& fine,
                    int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const T volfrac = T(1.0)/T(facx*facy);
    const int ii = i*facx;
    const int jj = j*facy;
    T c = 0;
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        c += fine(ii+iref,jj+jref,0,n+fcomp);
    }}
    crse(i,j,0,n+ccomp) = volfrac * c;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_nodes (Box const& bx, Array4<T> const& crse,
                          Array4<T const> const& fine,
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

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_nodes (int i, int j, int, int n, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    crse(i,j,0,n+ccomp) = fine(i*ratio[0],j*ratio[1],0,n+fcomp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_with_vol (int i, int j, int, int n, Array4<Real> const& crse,
                             Array4<Real const> const& fine,
                             Array4<Real const> const& fv,
                             int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int ii = i*facx;
    const int jj = j*facy;
    Real cd = 0.0, cv = 0.0;
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        cv +=                                 fv(ii+iref,jj+jref,0);
        cd += fine(ii+iref,jj+jref,0,n+fcomp)*fv(ii+iref,jj+jref,0);
    }}
    crse(i,j,0,n+ccomp) = cd/cv;
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

    for     (int n = 0; n < divu.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    divu(i,j,0,n) = dxi * (u(i+1,j,0,n)-u(i,j,0,n))
                        +           dyi * (v(i,j+1,0,n)-v(i,j,0,n));
                }
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

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_convective_difference (Box const& bx, Array4<amrex::Real> const& diff,
                                          Array4<Real const> const& u_face,
                                          Array4<Real const> const& v_face,
                                          Array4<Real const> const& s_on_x_face,
                                          Array4<Real const> const& s_on_y_face,
                                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for         (int n = 0; n <  diff.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    diff(i,j,n) = Real(0.5)*dxi * (u_face(i+1,j,0,0)+u_face(i,j,0,0)) *
                                            (s_on_x_face(i+1,j,0,n)-s_on_x_face(i,j,0,n))
                        +         Real(0.5)*dyi * (v_face(i,j+1,0,0)+v_face(i,j,0,0)) *
                                            (s_on_y_face(i,j+1,0,n)-s_on_y_face(i,j,0,n));
                }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence_rz (Box const& bx, Array4<Real> const& divu,
                                  Array4<Real const> const& u,
                                  Array4<Real const> const& v,
                                  Array4<Real const> const& ax,
                                  Array4<Real const> const& ay,
                                  Array4<Real const> const& vol) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for     (int n = 0; n < divu.ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    divu(i,j,0,n) = (ax(i+1,j,0,0)*u(i+1,j,0,n)-ax(i,j,0,0)*u(i,j,0,n))
                        +           (ay(i,j+1,0,0)*v(i,j+1,0,n)-ay(i,j,0,0)*v(i,j,0,n));
                    divu(i,j,0,n) /= vol(i,j,0,0);
                }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient_rz (Box const& bx, Array4<Real> const& grad,
                                Array4<Real const> const& u,
                                Array4<Real const> const& v,
                                Array4<Real const> const& ax,
                                Array4<Real const> const& ay,
                                Array4<Real const> const& vol) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            grad(i,j,0,0) = (ax(i+1,j,0,0)*u(i+1,j,0)-ax(i,j,0,0)*u(i,j,0))/vol(i,j,0,0);
            grad(i,j,0,1) = (ay(i,j+1,0,0)*v(i,j+1,0)-ay(i,j,0,0)*v(i,j,0))/vol(i,j,0,0);
        }
    }
}

}

#endif
