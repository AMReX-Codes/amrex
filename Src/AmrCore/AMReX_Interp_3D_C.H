#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST
inline
Vector<Real>
ccinterp_compute_voff (Box const& cbx, IntVect const& ratio, Geometry const& cgeom,
                       Geometry const& fgeom) noexcept
{
    const Box& fbx = amrex::refine(cbx,ratio);
    const auto& flen = amrex::length(fbx);
    const auto& flo = amrex::lbound(fbx);
    const auto& fhi = amrex::ubound(fbx);
    const int nxyz = flen.x + flen.y + flen.z;
    Vector<Real> voff(nxyz);
    
    const Real xrinv = 1./ratio[0];
    const Real yrinv = 1./ratio[1];
    const Real zrinv = 1./ratio[2];

    const Real xtmp = 0.5*(xrinv-1.0);
    const Real ytmp = 0.5*(yrinv-1.0);
    const Real ztmp = 0.5*(zrinv-1.0);

    Real* AMREX_RESTRICT xoff = &voff[0];
    AMREX_PRAGMA_SIMD
    for (int i = flo.x; i <= fhi.x; ++i) {
        const int ic = amrex::coarsen(i, ratio[0]);
        const int ioff = i - ic*ratio[0];
        xoff[i-flo.x] = ioff*xrinv + xtmp;
    }

    Real* AMREX_RESTRICT yoff = &voff[flen.x];
    AMREX_PRAGMA_SIMD
    for (int j = flo.y; j <= fhi.y; ++j) {
        const int jc = amrex::coarsen(j, ratio[1]);
        const int joff = j - jc*ratio[1];
        yoff[j-flo.y] = joff*yrinv + ytmp;
    }

    Real* AMREX_RESTRICT zoff = &voff[flen.x+flen.y];
    AMREX_PRAGMA_SIMD
    for (int k = flo.z; k <= fhi.z; ++k) {
        const int kc = amrex::coarsen(k, ratio[2]);
        const int koff = k - kc*ratio[2];
        zoff[k-flo.z] = koff*zrinv + ztmp;
    }

    return voff;
}

namespace {

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
compute_slopes (const Dim3& len, const Dim3& lo, const Dim3& hi, FabView<Real> const& slopes,
                FabView<const Real> const& u, Box const& sbox, const Dim3& slo, const Dim3& shi,
                const Dim3& slen, BCRec const& bc, const int n, const int ncomp) noexcept
{
    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,k,n        ) = 0.5*(u(i+1,j,k,n)-u(i-1,j,k,n));
                slopes(i,j,k,n+ncomp  ) = 0.5*(u(i,j+1,k,n)-u(i,j-1,k,n));
                slopes(i,j,k,n+ncomp*2) = 0.5*(u(i,j,k+1,n)-u(i,j,k-1,n));
            }
        }
    }

    if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
    {
        const int i = 0;
        if (slen.x >= 2) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    slopes(i,j,k,n) = -(16./15.)*u(i-1,j,k,n) + 0.5*u(i,j,k,n)
                        + (2./3.)*u(i+1,j,k,n) - 0.1*u(i+2,j,k,n);
                }
            }
        } else {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    slopes(i,j,k,n) = 0.25*(u(i+1,j,k,n)+5.*u(i,j,k,n)-6.*u(i-1,j,k,n));
                }
            }
        }
    }

    if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
    {
        const int i = len.x-1;
        if (slen.x >= 2) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    slopes(i,j,k,n) = (16./15.)*u(i+1,j,k,n) - 0.5*u(i,j,k,n)
                        - (2./3.)*u(i-1,j,k,n) + 0.1*u(i-2,j,k,n);
                }
            }                
        } else {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    slopes(i,j,k,n) = -0.25*(u(i-1,j,k,n)+5.*u(i,j,k,n)-6.*u(i+1,j,k,n));
                }
            }
        }
    }

    if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
    {
        const int j = 0;
        if (slen.y >= 2) {
            for     (int k = 0; k < len.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+ncomp) = -(16./15.)*u(i,j-1,k,n) + 0.5*u(i,j,k,n)
                        + (2./3.)*u(i,j+1,k,n) - 0.1*u(i,j+2,k,n);
                }
            }
        } else {
            for     (int k = 0; k < len.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+ncomp) = 0.25*(u(i,j+1,k,n)+5.*u(i,j,k,n)-6.*u(i,j-1,k,n));
                }
            }
        }
    }

    if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
    {
        const int j = len.y-1;
        if (slen.y >= 2) {
            for     (int k = 0; k < len.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+ncomp) = (16./15.)*u(i,j+1,k,n) - 0.5*u(i,j,k,n)
                        - (2./3.)*u(i,j-1,k,n) + 0.1*u(i,j-2,k,n);
                }
            }
        } else {
            for     (int k = 0; k < len.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+ncomp) = -0.25*(u(i,j-1,k,n)+5.*u(i,j,k,n)-6.*u(i,j+1,k,n));
                }
            }
        }
    }

    if (lo.z == slo.z && (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::hoextrap))
    {
        const int k = 0;
        if (slen.z >= 2) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+2*ncomp) = -(16./15.)*u(i,j,k-1,n) + 0.5*u(i,j,k,n)
                        + (2./3.)*u(i,j,k+1,n) - 0.1*u(i,j,k+2,n);
                }
            }
        } else {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+2*ncomp) = 0.25*(u(i,j,k+1,n)+5.*u(i,j,k,n)-6.*u(i,j,k-1,n));
                }
            }
        }
    }

    if (hi.z == shi.z && (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::hoextrap))
    {
        const int k = len.z-1;
        if (slen.z >= 2) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+2*ncomp) = (16./15.)*u(i,j,k+1,n) - 0.5*u(i,j,k,n)
                        - (2./3.)*u(i,j,k-1,n) + 0.1*u(i,j,k-2,n);
                }
            }
        } else {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n+2*ncomp) = -0.25*(u(i,j,k-1,n)+5.*u(i,j,k,n)-6.*u(i,j,k+1,n));
                }
            }
        }
    }
}

}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_linlim (Box const& bx, FArrayBox& ccfab,
                           FArrayBox const& ufab, const int icomp, const int ncomp,
                           BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto u = ufab.view(lo,icomp);

    const auto slopes = ccfab.view(lo,0);
    const auto sf     = ccfab.view(lo,ncomp*AMREX_SPACEDIM);  // slope factor

    const auto& sbox = ccfab.box();
    const auto slo  = amrex::lbound(sbox);
    const auto shi  = amrex::ubound(sbox);
    const auto slen = amrex::length(sbox);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                sf(i,j,k,0) = 1.0;
                sf(i,j,k,1) = 1.0;
                sf(i,j,k,2) = 1.0;
            }
        }
    }

    for (int n = 0; n < ncomp; ++n)
    {
        compute_slopes(len, lo, hi, slopes, u, sbox, slo, shi, slen, bcr[n], n, ncomp);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = 2.0*(u(i+1,j,k,n)-u(i  ,j,k,n));
                    Real back = 2.0*(u(i  ,j,k,n)-u(i-1,j,k,n));
                    Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,0) = amrex::min(sf(i,j,k,0), slopes(i,j,k,n)/cen);
                    } else {
                        sf(i,j,k,0) = 0.0;
                    }

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = 2.0*(u(i,j+1,k,n)-u(i,j  ,k,n));
                    back = 2.0*(u(i,j  ,k,n)-u(i,j-1,k,n));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,1) = amrex::min(sf(i,j,k,1), slopes(i,j,k,n+ncomp)/cen);
                    } else {
                        sf(i,j,k,1) = 0.0;
                    }

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = 2.0*(u(i,j,k+1,n)-u(i,j,k  ,n));
                    back = 2.0*(u(i,j,k  ,n)-u(i,j,k-1,n));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp*2) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,2) = amrex::min(sf(i,j,k,2), slopes(i,j,k,n+ncomp*2)/cen);
                    } else {
                        sf(i,j,k,2) = 0.0;
                    }
                }
            }
        }
    }

    for             (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slopes(i,j,k,n        ) *= sf(i,j,k,0);
                    slopes(i,j,k,n+ncomp  ) *= sf(i,j,k,1);
                    slopes(i,j,k,n+ncomp*2) *= sf(i,j,k,2);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_interp (Box const& bx, FArrayBox& finefab, const int fcomp, const int ncomp,
                    FArrayBox const& slopesfab, FArrayBox const& crsefab, const int ccomp,
                    Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto fine = finefab.view(lo,fcomp);

    const auto clo = amrex::lbound(slopesfab.box());
    const auto slopes = slopesfab.view(clo);
    const auto crse = crsefab.view(clo,ccomp);

    Box vbox = slopesfab.box();
    vbox.refine(ratio);
    const auto vlo  = amrex::lbound(vbox);
    const auto vlen = amrex::length(vbox);
    Real const* AMREX_RESTRICT xoff = &voff[lo.x-vlo.x];
    Real const* AMREX_RESTRICT yoff = &voff[lo.y-vlo.y+vlen.x];
    Real const* AMREX_RESTRICT zoff = &voff[lo.z-vlo.z+vlen.x+vlen.y];

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
            const int kc = amrex::coarsen(k+lo.z,ratio[2]) - clo.z;
            for (int j = 0; j < len.y; ++j) {
                const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                    fine(i,j,k,n) = crse(ic,jc,kc,n)
                        + xoff[i] * slopes(ic,jc,kc,n)
                        + yoff[j] * slopes(ic,jc,kc,n+ncomp)
                        + zoff[k] * slopes(ic,jc,kc,n+2*ncomp);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_mclim (Box const& bx, FArrayBox& ccfab,
                          FArrayBox const& ufab, const int icomp, const int ncomp,
                          BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto u = ufab.view(lo,icomp);

    const auto slopes = ccfab.view(lo,0);
    const auto mm     = ccfab.view(lo,ncomp*AMREX_SPACEDIM);  // min and max

    const auto& sbox = ccfab.box();
    const auto slo  = amrex::lbound(sbox);
    const auto shi  = amrex::ubound(sbox);
    const auto slen = amrex::length(sbox);

    for             (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    Real cmn = u(i,j,k,n);
                    Real cmx = cmn;
                    for         (int koff = -1; koff <= 1; ++koff) {
                        for     (int joff = -1; joff <= 1; ++joff) {
                            for (int ioff = -1; ioff <= 1; ++ioff) {
                                cmn = amrex::min(cmn,u(i+ioff,j+joff,k+koff,n));
                                cmx = amrex::max(cmx,u(i+ioff,j+joff,k+koff,n));
                            }
                        }
                    }
                    mm(i,j,k,n      ) = cmn - u(i,j,k,n);
                    mm(i,j,k,n+ncomp) = cmx - u(i,j,k,n);
                }
            }
        }

        compute_slopes(len, lo, hi, slopes, u, sbox, slo, shi, slen, bcr[n], n, ncomp);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = 2.0*(u(i+1,j,k,n)-u(i  ,j,k,n));
                    Real back = 2.0*(u(i  ,j,k,n)-u(i-1,j,k,n));
                    Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = 2.0*(u(i,j+1,k,n)-u(i,j  ,k,n));
                    back = 2.0*(u(i,j  ,k,n)-u(i,j-1,k,n));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = 2.0*(u(i,j,k+1,n)-u(i,j,k  ,n));
                    back = 2.0*(u(i,j,k  ,n)-u(i,j,k-1,n));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp*2) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_fine_alpha (Box const& bx, FArrayBox& alphafab, FArrayBox const& ccfab, const int ncomp,
                        Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto alpha = alphafab.view(lo);

    const auto clo = amrex::lbound(ccfab.box());
    const auto slopes = ccfab.view(clo);
    const auto mm     = ccfab.view(clo,ncomp*AMREX_SPACEDIM); // min and max

    const Box& alpha_box = alphafab.box();
    const auto vlo  = amrex::lbound(alpha_box);
    const auto vlen = amrex::length(alpha_box);
    Real const* AMREX_RESTRICT xoff = &voff[lo.x-vlo.x];
    Real const* AMREX_RESTRICT yoff = &voff[lo.y-vlo.y+vlen.x];
    Real const* AMREX_RESTRICT zoff = &voff[lo.z-vlo.z+vlen.x+vlen.y];

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
            const int kc = amrex::coarsen(k+lo.z,ratio[2]) - clo.z;
            for (int j = 0; j < len.y; ++j) {
                const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                    const Real dummy_fine = xoff[i]*slopes(ic,jc,kc,n)
                        +                   yoff[j]*slopes(ic,jc,kc,n+ncomp)
                        +                   zoff[k]*slopes(ic,jc,kc,n+ncomp*2);

                    if (dummy_fine > mm(ic,jc,kc,n+ncomp) && dummy_fine != 0.0) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n+ncomp) / dummy_fine;
                    } else if (dummy_fine < mm(ic,jc,kc,n) && dummy_fine != 0.0) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n) / dummy_fine;
                    } else {
                        alpha(i,j,k,n) = 1.0;
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_mmlim (Box const& bx, FArrayBox& ccfab, FArrayBox const& alphafab,
                          const int ncomp, IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto slopes = ccfab.view(lo);

    const auto flo = amrex::refine(lo,ratio);
    const auto alpha = alphafab.view(flo);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
            const int kk = k*ratio[2];
            for (int j = 0; j < len.y; ++j) {
                const int jj = j*ratio[1];
                for (int i = 0; i < len.x; ++i) {
                    const int ii = i*ratio[0];
                    Real a = 1.0;
                    for         (int koff = 0; koff < ratio[2]; ++koff) {
                        for     (int joff = 0; joff < ratio[1]; ++joff) {
                            for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                                a = amrex::min(a, alpha(ii+ioff,jj+joff,kk+koff,n));
                            }
                        }
                    }
                    slopes(i,j,k,n        ) *= a;
                    slopes(i,j,k,n+ncomp  ) *= a;
                    slopes(i,j,k,n+ncomp*2) *= a;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
pcinterp_interp (Box const& bx, FArrayBox& finefab, const int fcomp, const int ncomp,
                 FArrayBox const& crsefab, const int ccomp, IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto fine = finefab.view(lo,fcomp);

    const auto clo = amrex::lbound(crsefab.box());
    const auto crse = crsefab.view(clo,ccomp);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
            const int kc = amrex::coarsen(k+lo.z,ratio[2]) - clo.z;
            for (int j = 0; j < len.y; ++j) {
                const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                    fine(i,j,k,n) = crse(ic,jc,kc,n);
                }
            }
        }
    }
}

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int iz   = 2;
    static constexpr int ixy  = 3;
    static constexpr int ixz  = 4;
    static constexpr int iyz  = 5;
    static constexpr int ixyz = 6;
}

AMREX_GPU_HOST_DEVICE inline void
nodebilin_slopes (Box const& bx, FArrayBox& slopefab, FArrayBox const& ufab,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto slope = slopefab.view(lo);
    const auto u = ufab.view(lo,icomp);

    const Real rx = 1.0/ratio[0];
    const Real ry = 1.0/ratio[1];
    const Real rz = 1.0/ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    Real dx00 = u(i+1,j,k,n) - u(i,j,k,n);
                    Real d0x0 = u(i,j+1,k,n) - u(i,j,k,n);
                    Real d00x = u(i,j,k+1,n) - u(i,j,k,n);

                    Real dx10 = u(i+1,j+1,k,n) - u(i,j+1,k,n);
                    Real dx01 = u(i+1,j,k+1,n) - u(i,j,k+1,n);
                    Real d0x1 = u(i,j+1,k+1,n) - u(i,j,k+1,n);

                    Real dx11 = u(i+1,j+1,k+1,n) - u(i,j+1,k+1,n);

                    slope(i,j,k,n+ncomp*ix  ) = rx*dx00;
                    slope(i,j,k,n+ncomp*iy  ) = ry*d0x0;
                    slope(i,j,k,n+ncomp*iz  ) = rz*d00x;
                    slope(i,j,k,n+ncomp*ixy ) = rx*ry*(dx10 - dx00);
                    slope(i,j,k,n+ncomp*ixz ) = rx*rz*(dx01 - dx00);
                    slope(i,j,k,n+ncomp*iyz ) = ry*rz*(d0x1 - d0x0);
                    slope(i,j,k,n+ncomp*ixyz) = rx*ry*rz*(dx11 - dx01 - dx10 + dx00);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
nodebilin_interp (Box const& bx, FArrayBox& finefab, const int fcomp, const int ncomp,
                  FArrayBox const& slopefab, FArrayBox const& crsefab, const int ccomp,
                  IntVect const& ratio) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto fine = finefab.view(lo,fcomp);

    const Box& sbox = slopefab.box();
    const auto clo = amrex::lbound(sbox);
    const auto chi = amrex::ubound(sbox);
    const auto slope = slopefab.view(clo);
    const auto crse = crsefab.view(clo,ccomp);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
            const int kc = amrex::min(amrex::coarsen(k+lo.z,ratio[2]),chi.z) - clo.z;
            const Real fz = (k+lo.z) - (kc+clo.z)*ratio[2];
            for (int j = 0; j < len.y; ++j) {
                const int jc = amrex::min(amrex::coarsen(j+lo.y,ratio[1]),chi.y) - clo.y;
                const Real fy = (j+lo.y) - (jc+clo.y)*ratio[1];
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    const int ic = amrex::min(amrex::coarsen(i+lo.x,ratio[0]),chi.x) - clo.x;
                    const Real fx = (i+lo.x) - (ic+clo.x)*ratio[0];
                    fine(i,j,k,n) = crse(ic,jc,kc,n)
                        + fx*slope(ic,jc,kc,n+ncomp*ix)
                        + fy*slope(ic,jc,kc,n+ncomp*iy)
                        + fz*slope(ic,jc,kc,n+ncomp*iz)
                        + fx*fy*slope(ic,jc,kc,n+ncomp*ixy)
                        + fx*fz*slope(ic,jc,kc,n+ncomp*ixz)
                        + fy*fz*slope(ic,jc,kc,n+ncomp*iyz)
                        + fx*fy*fz*slope(ic,jc,kc,n+ncomp*ixyz);
                }
            }
        }
    }
}

}

#endif
