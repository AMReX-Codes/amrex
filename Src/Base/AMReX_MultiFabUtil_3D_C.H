#ifndef AMREX_MULTIFAB_UTIL_3D_C_H_
#define AMREX_MULTIFAB_UTIL_3D_C_H_

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
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            cc(i,j,k,n+cccomp) = 0.125*( nd(i,j  ,k  ,n+ndcomp) + nd(i+1,j  ,k  ,n+ndcomp)
                                       + nd(i,j+1,k  ,n+ndcomp) + nd(i+1,j+1,k  ,n+ndcomp)
                                       + nd(i,j  ,k+1,n+ndcomp) + nd(i+1,j  ,k+1,n+ndcomp)
                                       + nd(i,j+1,k+1,n+ndcomp) + nd(i+1,j+1,k+1,n+ndcomp));
        }}}
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_eg_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& Ex,
                         Array4<Real const> const& Ey,
                         Array4<Real const> const& Ez,
                         int cccomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        cc(i,j,k,0+cccomp) = 0.25 * ( Ex(i,j,k) + Ex(i,j+1,k) + Ex(i,j,k+1) + Ex(i,j+1,k+1) );
        cc(i,j,k,1+cccomp) = 0.25 * ( Ey(i,j,k) + Ey(i+1,j,k) + Ey(i,j,k+1) + Ey(i+1,j,k+1) );
        cc(i,j,k,2+cccomp) = 0.25 * ( Ez(i,j,k) + Ez(i+1,j,k) + Ez(i,j+1,k) + Ez(i+1,j+1,k) );
    }}}
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_fc_to_cc (Box const& bx,
                         Array4<Real      > const& cc,
                         Array4<Real const> const& fx,
                         Array4<Real const> const& fy,
                         Array4<Real const> const& fz,
                         int cccomp) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        cc(i,j,k,0+cccomp) = 0.5 * ( fx(i,j,k) + fx(i+1,j,k) );
        cc(i,j,k,1+cccomp) = 0.5 * ( fy(i,j,k) + fy(i,j+1,k) );
        cc(i,j,k,2+cccomp) = 0.5 * ( fz(i,j,k) + fz(i,j,k+1) );
    }}}
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, Box const& ybx, Box const& zbx,
                         Array4<Real> const& fx, Array4<Real> const& fy,
                         Array4<Real> const& fz, Array4<Real const> const& cc) noexcept
{
    const auto ndlo = lbound(ndbx);
    const auto ndhi = ubound(ndbx);
    const auto xlo = amrex::elemwiseMax(ndlo, lbound(xbx));
    const auto ylo = amrex::elemwiseMax(ndlo, lbound(ybx));
    const auto zlo = amrex::elemwiseMax(ndlo, lbound(zbx));
    const auto xhi = amrex::elemwiseMin(ndhi, ubound(xbx));
    const auto yhi = amrex::elemwiseMin(ndhi, ubound(ybx));
    const auto zhi = amrex::elemwiseMin(ndhi, ubound(zbx));

    for         (int k = xlo.z; k <= xhi.z; ++k) {
        for     (int j = xlo.y; j <= xhi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = xlo.x; i <= xhi.x; ++i) {
                fx(i,j,k) = 0.5*(cc(i-1,j,k) + cc(i,j,k));
            }
        }
    }    

    for         (int k = ylo.z; k <= yhi.z; ++k) {
        for     (int j = ylo.y; j <= yhi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = ylo.x; i <= yhi.x; ++i) {
                fy(i,j,k) = 0.5*(cc(i,j-1,k) + cc(i,j,k));
            }
        }
    }    

    for         (int k = zlo.z; k <= zhi.z; ++k) {
        for     (int j = zlo.y; j <= zhi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = zlo.x; i <= zhi.x; ++i) {
                fz(i,j,k) = 0.5*(cc(i,j,k-1) + cc(i,j,k));
            }
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
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / (facy*facz);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
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
        Real facInv = 1.0 / (facx*facz);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
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
        Real facInv = 1.0 / (facx*facy);
        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
            for (int j = clo.y; j <= chi.y; ++j) {
            for (int i = clo.x; i <= chi.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                int kk = k*facz;
                Real c = 0.;
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
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facx;
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
        Real facInv = 1.0 / facy;
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
        Real facInv = 1.0 / facz;
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
    const int facz = ratio[2];
    const Real volfrac = 1.0/static_cast<Real>(facx*facy*facz);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = clo.z; k <= chi.z; ++k) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            int kk = k*facz;
            Real c = 0.;
            for (int kref = 0; kref < facz; ++kref) {
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                c += fine(ii+iref,jj+jref,kk+kref,n+fcomp);
            }}}
            crse(i,j,k,n+ccomp) = volfrac * c;
        }}}
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
    const int facz = ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        for (int k = clo.z; k <= chi.z; ++k) {
        for (int j = clo.y; j <= chi.y; ++j) {
        for (int i = clo.x; i <= chi.x; ++i) {
            int ii = i*facx;
            int jj = j*facy;
            int kk = k*facz;
            Real cd = 0., cv = 0.;
            for (int kref = 0; kref < facz; ++kref) {
            for (int jref = 0; jref < facy; ++jref) {
            for (int iref = 0; iref < facx; ++iref) {
                cv +=                                       fv(ii+iref,jj+jref,kk+kref);
                cd += fine(ii+iref,jj+jref,kk+kref,n+fcomp)*fv(ii+iref,jj+jref,kk+kref);
            }}}
            crse(i,j,k,n+ccomp) = cd/cv;
        }}}
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

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                divu(i,j,k) = dxi * (u(i+1,j,k)-u(i,j,k))
                    +         dyi * (v(i,j+1,k)-v(i,j,k))
                    +         dzi * (w(i,j,k+1)-w(i,j,k));
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

} // namespace amrex

#endif
