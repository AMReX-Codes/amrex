#ifndef AMREX_MULTIFAB_UTIL_2D_C_H_
#define AMREX_MULTIFAB_UTIL_2D_C_H_

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
        for (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,0,n) = 0.25*(nd(i,j,0,n)+nd(i+1,j,0,n)+nd(i,j+1,0,n)+nd(i+1,j+1,0,n));
        }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& exfab, FArrayBox const& eyfab,
                         int cccomp) noexcept
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
inline
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         FArrayBox const& fxfab, FArrayBox const& fyfab,
                         int cccomp) noexcept
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
inline
void amrex_avg_cc_to_fc (Box const& ndbx, Box const& xbx, Box const& ybx,
                         FArrayBox& fxfab, FArrayBox& fyfab,
                         FArrayBox const& ccfab) noexcept
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

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facy;
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    Real c = 0.;
                    for (int jref = 0; jref < facy; ++jref) {
                        c += fine(ii,jj+jref,0,n);
                    }
                    crse(i,j,0,n) = c * facInv;
                }
            }
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / facx;
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    Real c = 0.;
                    for (int iref = 0; iref < facx; ++iref) {
                        c += fine(ii+iref,jj,0,n);
                    }
                    crse(i,j,0,n) = c * facInv;
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

    switch (idir) {
    case 0:
    {
        Real facInv = 1.0 / facx;
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    Real c = 0.;
                    for (int iref = 0; iref < facx; ++iref) {
                        c += fine(ii+iref,jj,0,n);
                    }
                    crse(i,j,0,n) = c * facInv;
                }
            }
        }
        break;
    }
    case 1:
    {
        Real facInv = 1.0 / facy;
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int ii = i*facx;
                    int jj = j*facy;
                    Real c = 0.;
                    for (int jref = 0; jref < facx; ++jref) {
                        c += fine(ii,jj+jref,0,n);
                    }
                    crse(i,j,0,n) = c * facInv;
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
    const Real volfrac = 1.0/static_cast<Real>(facx*facy);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real c = 0.;
                for     (int jref = 0; jref < facy; ++jref) {
                    for (int iref = 0; iref < facx; ++iref) {
                        c += fine(ii+iref,jj+jref,0,n);
                    }
                }
                crse(i,j,0,n) = volfrac * c;
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

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            int jj = j*facy;
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                crse(i,j,0,n) = fine(i*facx,jj,0,n);
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

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                int ii = i*facx;
                int jj = j*facy;
                Real cd = 0., cv = 0.;
                for     (int jref = 0; jref < facy; ++jref) {
                    for (int iref = 0; iref < facx; ++iref) {
                        cv +=                           fv(ii+iref,jj+jref,0);
                        cd += fine(ii+iref,jj+jref,0,n)*fv(ii+iref,jj+jref,0);
                    }
                }
                crse(i,j,0,n) = cd/cv;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void amrex_compute_divergence (Box const& bx, FArrayBox& divufab,
                               FArrayBox const& ufab, FArrayBox const& vfab,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto divu = divufab.view(lo);
    const auto    u =    ufab.view(lo);
    const auto    v =    vfab.view(lo);

    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];

    for     (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
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
