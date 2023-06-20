#ifndef AMREX_MULTIFAB_UTIL_3D_C_H_
#define AMREX_MULTIFAB_UTIL_3D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_nd_to_cc (int i, int j, int k, int n,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& nd,
                         int cccomp, int ndcomp) noexcept
{
    cc(i,j,k,n+cccomp) = Real(0.125)*( nd(i,j  ,k  ,n+ndcomp) + nd(i+1,j  ,k  ,n+ndcomp)
                                     + nd(i,j+1,k  ,n+ndcomp) + nd(i+1,j+1,k  ,n+ndcomp)
                                     + nd(i,j  ,k+1,n+ndcomp) + nd(i+1,j  ,k+1,n+ndcomp)
                                     + nd(i,j+1,k+1,n+ndcomp) + nd(i+1,j+1,k+1,n+ndcomp));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_eg_to_cc (int i, int j, int k,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         Array4<Real const> const& Ey,
                         Array4<Real const> const& Ez,
                         int cccomp) noexcept
{
    cc(i,j,k,0+cccomp) = Real(0.25) * ( Ex(i,j,k) + Ex(i,j+1,k) + Ex(i,j,k+1) + Ex(i,j+1,k+1) );
    cc(i,j,k,1+cccomp) = Real(0.25) * ( Ey(i,j,k) + Ey(i+1,j,k) + Ey(i,j,k+1) + Ey(i+1,j,k+1) );
    cc(i,j,k,2+cccomp) = Real(0.25) * ( Ez(i,j,k) + Ez(i+1,j,k) + Ez(i,j+1,k) + Ez(i+1,j+1,k) );
}

template <typename CT, typename FT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_fc_to_cc (int i, int j, int k,
                         Array4<CT      > const& cc,
                         Array4<FT const> const& fx,
                         Array4<FT const> const& fy,
                         Array4<FT const> const& fz,
                         int cccomp) noexcept
{
    cc(i,j,k,0+cccomp) = CT(0.5) * CT( fx(i,j,k) + fx(i+1,j,k) );
    cc(i,j,k,1+cccomp) = CT(0.5) * CT( fy(i,j,k) + fy(i,j+1,k) );
    cc(i,j,k,2+cccomp) = CT(0.5) * CT( fz(i,j,k) + fz(i,j,k+1) );
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avg_cc_to_fc (int i, int j, int k, int n, Box const& xbx, Box const& ybx, Box const& zbx,
                         Array4<Real> const& fx, Array4<Real> const& fy,
                         Array4<Real> const& fz, Array4<Real const> const& cc,
                         bool use_harmonic_averaging) noexcept
{
    if (use_harmonic_averaging)
    {
        if (xbx.contains(i,j,k)) {
            if (cc(i-1,j,k,n) == Real(0.0) || cc(i,j,k,n) == Real(0.0)) {
                fx(i,j,k,n) = Real(0.0);
            } else {
                fx(i,j,k,n) = Real(2.0) / (Real(1.0)/cc(i-1,j,k,n) + Real(1.0)/cc(i,j,k,n));
            }
        }
        if (ybx.contains(i,j,k)) {
            if (cc(i,j-1,k,n) == Real(0.0) || cc(i,j,k,n) == Real(0.0)) {
                fy(i,j,k,n) = Real(0.0);
            } else {
                fy(i,j,k,n) = Real(2.0) / (Real(1.0)/cc(i,j-1,k,n) + Real(1.0)/cc(i,j,k,n));
            }
        }
        if (zbx.contains(i,j,k)) {
            if (cc(i,j,k-1,n) == Real(0.0) || cc(i,j,k,n) == Real(0.0)) {
                fz(i,j,k,n) = Real(0.0);
            } else {
                fz(i,j,k,n) = Real(2.0) / (Real(1.0)/cc(i,j,k-1,n) + Real(1.0)/cc(i,j,k,n));
            }
        }
    } else {
        if (xbx.contains(i,j,k)) {
            fx(i,j,k,n) = Real(0.5)*(cc(i-1,j,k,n) + cc(i,j,k,n));
        }
        if (ybx.contains(i,j,k)) {
            fy(i,j,k,n) = Real(0.5)*(cc(i,j-1,k,n) + cc(i,j,k,n));
        }
        if (zbx.contains(i,j,k)) {
            fz(i,j,k,n) = Real(0.5)*(cc(i,j,k-1,n) + cc(i,j,k,n));
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
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        T facInv = T(1.0) / (facy*facz);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                T c = T(0.);
                for (int kref = 0; kref < facz; ++kref) {
                for (int jref = 0; jref < facy; ++jref) {
                    c += fine(ii,jj+jref,kk+kref,n+fcomp);
                }}
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    case 1:
    {
        T facInv = T(1.0) / (facx*facz);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                T c = T(0.);
                for (int kref = 0; kref < facz; ++kref) {
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,kk+kref,n+fcomp);
                }}
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    case 2:
    {
        T facInv = T(1.0) / (facx*facy);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                T c = T(0.);
                for (int jref = 0; jref < facy; ++jref) {
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj+jref,kk,n+fcomp);
                }}
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_faces (int i, int j, int k, int n, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;

    switch (idir) {
    case 0:
    {
        const T facInv = T(1.0) / (facy*facz);
        T c = T(0.);
        for (int kref = 0; kref < facz; ++kref) {
        for (int jref = 0; jref < facy; ++jref) {
            c += fine(ii,jj+jref,kk+kref,n+fcomp);
        }}
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const T facInv = T(1.0) / (facx*facz);
        T c = T(0.);
        for (int kref = 0; kref < facz; ++kref) {
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,kk+kref,n+fcomp);
        }}
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
    case 2:
    {
        const T facInv = T(1.0) / (facx*facy);
        T c = T(0.);
        for (int jref = 0; jref < facy; ++jref) {
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj+jref,kk,n+fcomp);
        }}
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
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
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        Real facInv = Real(1.0) / static_cast<Real>(facx);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
                for (int iref = 0; iref < facx; ++iref) {
                    c += fine(ii+iref,jj,kk,n+fcomp);
                }
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    case 1:
    {
        Real facInv = Real(1.0) / static_cast<Real>(facy);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
                for (int jref = 0; jref < facy; ++jref) {
                    c += fine(ii,jj+jref,kk,n+fcomp);
                }
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    case 2:
    {
        Real facInv = Real(1.0) / static_cast<Real>(facz);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
                for (int kref = 0; kref < facz; ++kref) {
                    c += fine(ii,jj,kk+kref,n+fcomp);
                }
                crse(i,j,k,n+ccomp) = c * facInv;
            }}}
        }
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_edges (int i, int j, int k, int n, Array4<Real> const& crse,
                          Array4<Real const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio, int idir) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;

    switch (idir) {
    case 0:
    {
        const Real facInv = Real(1.0) / static_cast<Real>(facx);
        Real c = 0.;
        for (int iref = 0; iref < facx; ++iref) {
            c += fine(ii+iref,jj,kk,n+fcomp);
        }
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
    case 1:
    {
        const Real facInv = Real(1.0) / static_cast<Real>(facy);
        Real c = 0.;
        for (int jref = 0; jref < facy; ++jref) {
            c += fine(ii,jj+jref,kk,n+fcomp);
        }
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
    case 2:
    {
        const Real facInv = Real(1.0) / static_cast<Real>(facz);
        Real c = 0.;
        for (int kref = 0; kref < facz; ++kref) {
            c += fine(ii,jj,kk+kref,n+fcomp);
        }
        crse(i,j,k,n+ccomp) = c * facInv;
        break;
    }
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
    const int facz = ratio[2];
    const T volfrac = T(1.0)/T(facx*facy*facz);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = clo.z; k <= chi.z; ++k) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            int kk = k*facz;
            T c = 0;
            for (int kref = 0; kref < facz; ++kref) {
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,kk+kref,n+fcomp);
            }}}
            crse(i,j,k,n+ccomp) = volfrac * c;
        }}}
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown (int i, int j, int k, int n, Array4<T> const& crse,
                    Array4<T const> const& fine,
                    int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const T volfrac = T(1.0)/T(facx*facy*facz);
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;
    T c = 0;
    for (int kref = 0; kref < facz; ++kref) {
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        c += fine(ii+iref,jj+jref,kk+kref,n+fcomp);
    }}}
    crse(i,j,k,n+ccomp) = volfrac * c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_with_vol (int i, int j, int k, int n, Array4<Real> const& crse,
                             Array4<Real const> const& fine,
                             Array4<Real const> const& fv,
                             int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;
    Real cd = 0.0, cv = 0.0;
    for (int kref = 0; kref < facz; ++kref) {
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        cv +=                                       fv(ii+iref,jj+jref,kk+kref);
        cd += fine(ii+iref,jj+jref,kk+kref,n+fcomp)*fv(ii+iref,jj+jref,kk+kref);
    }}}
    crse(i,j,k,n+ccomp) = cd/cv;
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
    const int facz = ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = clo.z; k <= chi.z; ++k) {
            int kk = k*facz;
            for     (int j = clo.y; j <= chi.y; ++j) {
                int jj = j*facy;
                AMREX_PRAGMA_SIMD
                for (int i = clo.x; i <= chi.x; ++i) {
                    crse(i,j,k,n+ccomp) = fine(i*facx,jj,kk,n+fcomp);
                }
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_avgdown_nodes (int i, int j, int k, int n, Array4<T> const& crse,
                          Array4<T const> const& fine,
                          int ccomp, int fcomp, IntVect const& ratio) noexcept
{
    crse(i,j,k,n+ccomp) = fine(i*ratio[0],j*ratio[1],k*ratio[2],n+fcomp);
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, Array4<Real> const& divu,
                               Array4<Real const> const& u,
                               Array4<Real const> const& v,
                               Array4<Real const> const& w,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];

    for             (int n = 0; n <  divu.ncomp; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        divu(i,j,k,n) = dxi * (u(i+1,j,k,n)-u(i,j,k,n))
                            +           dyi * (v(i,j+1,k,n)-v(i,j,k,n))
                            +           dzi * (w(i,j,k+1,n)-w(i,j,k,n));
                    }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_gradient (Box const& bx, Array4<Real> const& grad,
                             Array4<Real const> const& u,
                             Array4<Real const> const& v,
                             Array4<Real const> const& w,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                grad(i,j,k,0) = dxi * (u(i+1,j,k,0)-u(i,j,k,0));
                grad(i,j,k,1) = dyi * (v(i,j+1,k,0)-v(i,j,k,0));
                grad(i,j,k,2) = dzi * (w(i,j,k+1,0)-w(i,j,k,0));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_convective_difference (Box const& bx, Array4<Real> const& diff,
                                          Array4<Real const> const& u_face,
                                          Array4<Real const> const& v_face,
                                          Array4<Real const> const& w_face,
                                          Array4<Real const> const& s_on_x_face,
                                          Array4<Real const> const& s_on_y_face,
                                          Array4<Real const> const& s_on_z_face,
                                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];

    for             (int n = 0; n <  diff.ncomp; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        diff(i,j,k,n) = Real(0.5)*dxi * (u_face(i+1,j,k)+u_face(i,j,k)) *
                                                  (s_on_x_face(i+1,j,k,n)-s_on_x_face(i,j,k,n))
                            +           Real(0.5)*dyi * (v_face(i,j+1,k)+v_face(i,j,k)) *
                                                  (s_on_y_face(i,j+1,k,n)-s_on_y_face(i,j,k,n))
                            +           Real(0.5)*dzi * (w_face(i,j,k+1)+w_face(i,j,k)) *
                                                  (s_on_z_face(i,j,k+1,n)-s_on_z_face(i,j,k,n));
                    }
            }
        }
    }
}

} // namespace amrex

#endif
