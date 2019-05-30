#ifndef AMREX_MULTIFAB_UTIL_3D_C_H_
#define AMREX_MULTIFAB_UTIL_3D_C_H_

#include <AMReX_Gpu.H>
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
inline
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& exfab, FArrayBox const& eyfab, FArrayBox const& ezfab,
                         int cccomp) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
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
inline
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, FArrayBox const& fyfab, FArrayBox const& fzfab,
                         int cccomp) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
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

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, Box const& ybx, Box const& zbx,
                         FArrayBox& fxfab, FArrayBox& fyfab, FArrayBox& fzfab,
                         FArrayBox const& ccfab) noexcept
{
    const auto lo = lbound(ndbx);
    const auto fx = fxfab.view(lo);
    const auto fy = fyfab.view(lo);
    const auto fz = fzfab.view(lo);
    const auto cc = ccfab.view(lo);

    const auto xlen = length(ndbx,xbx);
    for         (int k = 0; k < xlen.z; ++k) {
        for     (int j = 0; j < xlen.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < xlen.x; ++i) {
                fx(i,j,k) = 0.5*(cc(i-1,j,k) + cc(i,j,k));
            }
        }
    }    

    const auto ylen = length(ndbx,ybx);
    for         (int k = 0; k < ylen.z; ++k) {
        for     (int j = 0; j < ylen.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < ylen.x; ++i) {
                fy(i,j,k) = 0.5*(cc(i,j-1,k) + cc(i,j,k));
            }
        }
    }    

    const auto zlen = length(ndbx,zbx);
    for         (int k = 0; k < zlen.z; ++k) {
        for     (int j = 0; j < zlen.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < zlen.x; ++i) {
                fz(i,j,k) = 0.5*(cc(i,j,k-1) + cc(i,j,k));
            }
        }
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
    const int facy = ratio[1];
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / (facy*facz);
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for     (int kref = 0; kref < facz; ++kref) {
                            for (int jref = 0; jref < facy; ++jref) {
                                c += fine(ii,jj+jref,kk+kref,n);
                            }
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / (facx*facz);
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for     (int kref = 0; kref < facz; ++kref) {
                            for (int iref = 0; iref < facx; ++iref) {
                                c += fine(ii+iref,jj,kk+kref,n);
                            }
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
    }
    case 2:
    {
        Real facInv = 1.0 / (facx*facy);
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for     (int jref = 0; jref < facy; ++jref) {
                            for (int iref = 0; iref < facx; ++iref) {
                                c += fine(ii+iref,jj+jref,kk,n);
                            }
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
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
    const int facy = ratio[1];
    const int facz = ratio[2];

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facx;
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for (int iref = 0; iref < facx; ++iref) {
                            c += fine(ii+iref,jj,kk,n);
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / facy;
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for (int jref = 0; jref < facy; ++jref) {
                            c += fine(ii,jj+jref,kk,n);
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
    }
    case 2:
    {
        Real facInv = 1.0 / facz;
        for (int n = 0; n < ncomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        int ii = i*facx;
                        int jj = j*facy;
                        int kk = k*facz;
                        Real c = 0.;
                        for (int kref = 0; kref < facz; ++kref) {
                            c += fine(ii,jj,kk+kref,n);
                        }
                        crse(i,j,k,n) = c * facInv;
                    }
                }
            }
        }
        break;
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
    const int facy = ratio[1];
    const int facz = ratio[2];
    const Real volfrac = 1.0/static_cast<Real>(facx*facy*facz);

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    int kk = k*facz;
                    Real c = 0.;
                    for         (int kref = 0; kref < facz; ++kref) {
                        for     (int jref = 0; jref < facy; ++jref) {
                            for (int iref = 0; iref < facx; ++iref) {
                                c += fine(ii+iref,jj+jref,kk+kref,n);
                            }
                        }
                    }
                    crse(i,j,k,n) = volfrac * c;
                }
            }
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
    const int facy = ratio[1];
    const int facz = ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    int kk = k*facz;
                    Real cd = 0., cv = 0.;
                    for         (int kref = 0; kref < facz; ++kref) {
                        for     (int jref = 0; jref < facy; ++jref) {
                            for (int iref = 0; iref < facx; ++iref) {
                                cv +=                                 fv(ii+iref,jj+jref,kk+kref);
                                cd += fine(ii+iref,jj+jref,kk+kref,n)*fv(ii+iref,jj+jref,kk+kref);
                            }
                        }
                    }
                    crse(i,j,k,n) = cd/cv;
                }
            }
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
    const int facy = ratio[1];
    const int facz = ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            int kk = k*facz;
            for     (int j = 0; j < len.y; ++j) {
                int jj = j*facy;
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    crse(i,j,k,n) = fine(i*facx,jj,kk,n);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, FArrayBox& divufab,
                               FArrayBox const& ufab, FArrayBox const& vfab, FArrayBox const& wfab,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto divu = divufab.view(lo);
    const auto    u =    ufab.view(lo);
    const auto    v =    vfab.view(lo);
    const auto    w =    wfab.view(lo);

    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
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
