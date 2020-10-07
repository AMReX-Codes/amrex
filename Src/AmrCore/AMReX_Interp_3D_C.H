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
ccinterp_compute_voff (Box const& cbx, IntVect const& ratio, Geometry const& /*cgeom*/,
                       Geometry const& /*fgeom*/) noexcept
{
    const Box& fbx = amrex::refine(cbx,ratio);
    const auto& flen = amrex::length(fbx);
    const auto& flo = amrex::lbound(fbx);
    const auto& fhi = amrex::ubound(fbx);
    const int nxyz = flen.x + flen.y + flen.z;
    Vector<Real> voff(nxyz);
    
    const Real xrinv = Real(1.)/ratio[0];
    const Real yrinv = Real(1.)/ratio[1];
    const Real zrinv = Real(1.)/ratio[2];

    const Real xtmp = Real(0.5)*(xrinv-Real(1.));
    const Real ytmp = Real(0.5)*(yrinv-Real(1.));
    const Real ztmp = Real(0.5)*(zrinv-Real(1.));

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
compute_slopes (const Dim3& lo, const Dim3& hi,
                Array4<Real> const& slopes, const Dim3& slo, const Dim3& shi, int ns,
                Array4<Real const> const& u, int nu,
                BCRec const& bc, int ncomp) noexcept
{
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,k,ns        ) = Real(0.5)*(u(i+1,j,k,nu)-u(i-1,j,k,nu));
                slopes(i,j,k,ns+ncomp  ) = Real(0.5)*(u(i,j+1,k,nu)-u(i,j-1,k,nu));
                slopes(i,j,k,ns+ncomp*2) = Real(0.5)*(u(i,j,k+1,nu)-u(i,j,k-1,nu));
            }
        }
    }

    if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
    {
        const int i = slo.x;
        if (shi.x-slo.x >= 1) {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = -Real(16./15.)*u(i-1,j,k,nu) + Real(0.5)*u(i,j,k,nu)
                        + Real(2./3.)*u(i+1,j,k,nu) - Real(0.1)*u(i+2,j,k,nu);
                }
            }
        } else {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = Real(0.25)*(u(i+1,j,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i-1,j,k,nu));
                }
            }
        }
    }

    if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
    {
        const int i = shi.x;
        if (shi.x-slo.x >= 1) {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = Real(16./15.)*u(i+1,j,k,nu) - Real(0.5)*u(i,j,k,nu)
                        - Real(2./3.)*u(i-1,j,k,nu) + Real(0.1)*u(i-2,j,k,nu);
                }
            }                
        } else {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = -Real(0.25)*(u(i-1,j,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i+1,j,k,nu));
                }
            }
        }
    }

    if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
    {
        const int j = slo.y;
        if (shi.y-slo.y >= 1) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = -Real(16./15.)*u(i,j-1,k,nu) + Real(0.5)*u(i,j,k,nu)
                        + Real(2./3.)*u(i,j+1,k,nu) - Real(0.1)*u(i,j+2,k,nu);
                }
            }
        } else {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = Real(0.25)*(u(i,j+1,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j-1,k,nu));
                }
            }
        }
    }

    if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
    {
        const int j = shi.y;
        if (shi.y-slo.y >= 1) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = Real(16./15.)*u(i,j+1,k,nu) - Real(0.5)*u(i,j,k,nu)
                        - Real(2./3.)*u(i,j-1,k,nu) + Real(0.1)*u(i,j-2,k,nu);
                }
            }
        } else {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = -Real(0.25)*(u(i,j-1,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j+1,k,nu));
                }
            }
        }
    }

    if (lo.z == slo.z && (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::hoextrap))
    {
        const int k = slo.z;
        if (shi.z-slo.z >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = -Real(16./15.)*u(i,j,k-1,nu) + Real(0.5)*u(i,j,k,nu)
                        + Real(2./3.)*u(i,j,k+1,nu) - Real(0.1)*u(i,j,k+2,nu);
                }
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = Real(0.25)*(u(i,j,k+1,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j,k-1,nu));
                }
            }
        }
    }

    if (hi.z == shi.z && (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::hoextrap))
    {
        const int k = shi.z;
        if (shi.z-slo.z >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = Real(16./15.)*u(i,j,k+1,nu) - Real(0.5)*u(i,j,k,nu)
                        - Real(2./3.)*u(i,j,k-1,nu) + Real(0.1)*u(i,j,k-2,nu);
                }
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = -Real(0.25)*(u(i,j,k-1,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j,k+1,nu));
                }
            }
        }
    }
}

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cellconslin_slopes_linlim (Box const& bx, Array4<Real> const& slopes,
                           Array4<Real const> const& u, const int icomp, const int ncomp,
                           BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const auto slo  = amrex::lbound(slopes);
    const auto shi  = amrex::ubound(slopes);

    const Array4<Real> sf(slopes, ncomp*AMREX_SPACEDIM);  // slope factor

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                sf(i,j,k,0) = Real(1.);
                sf(i,j,k,1) = Real(1.);
                sf(i,j,k,2) = Real(1.);
            }
        }
    }

    for (int n = 0; n < ncomp; ++n)
    {
        const int nu = n + icomp;
        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = Real(2.)*(u(i+1,j,k,nu)-u(i  ,j,k,nu));
                    Real back = Real(2.)*(u(i  ,j,k,nu)-u(i-1,j,k,nu));
                    Real slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                    if (cen != Real(0.)) {
                        sf(i,j,k,0) = amrex::min(sf(i,j,k,0), slopes(i,j,k,n)/cen);
                    } else {
                        sf(i,j,k,0) = Real(0.);
                    }

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = Real(2.)*(u(i,j+1,k,nu)-u(i,j  ,k,nu));
                    back = Real(2.)*(u(i,j  ,k,nu)-u(i,j-1,k,nu));
                    slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n+ncomp) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                    if (cen != Real(0.)) {
                        sf(i,j,k,1) = amrex::min(sf(i,j,k,1), slopes(i,j,k,n+ncomp)/cen);
                    } else {
                        sf(i,j,k,1) = Real(0.);
                    }

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = Real(2.)*(u(i,j,k+1,nu)-u(i,j,k  ,nu));
                    back = Real(2.)*(u(i,j,k  ,nu)-u(i,j,k-1,nu));
                    slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n+ncomp*2) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                    if (cen != Real(0.)) {
                        sf(i,j,k,2) = amrex::min(sf(i,j,k,2), slopes(i,j,k,n+ncomp*2)/cen);
                    } else {
                        sf(i,j,k,2) = Real(0.);
                    }
                }
            }
        }
    }

    for             (int n = 0; n < ncomp; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,n        ) *= sf(i,j,k,0);
                    slopes(i,j,k,n+ncomp  ) *= sf(i,j,k,1);
                    slopes(i,j,k,n+ncomp*2) *= sf(i,j,k,2);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cellconslin_interp (Box const& bx,
                    Array4<Real> const& fine, const int fcomp, const int ncomp,
                    Array4<Real const> const& slopes,
                    Array4<Real const> const& crse, const int ccomp,
                    Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    Box vbox(slopes);
    vbox.refine(ratio);
    const auto vlo  = amrex::lbound(vbox);
    const auto vlen = amrex::length(vbox);
    Real const* AMREX_RESTRICT xoff = voff;
    Real const* AMREX_RESTRICT yoff = voff + vlen.x;
    Real const* AMREX_RESTRICT zoff = voff + (vlen.x+vlen.y);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp)
                        + xoff[i-vlo.x] * slopes(ic,jc,kc,n)
                        + yoff[j-vlo.y] * slopes(ic,jc,kc,n+ncomp)
                        + zoff[k-vlo.z] * slopes(ic,jc,kc,n+2*ncomp);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cellconslin_slopes_mclim (Box const& bx, Array4<Real> const& slopes,
                          Array4<Real const> const& u, const int icomp, const int ncomp,
                          BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const auto slo = amrex::lbound(slopes);
    const auto shi = amrex::ubound(slopes);

    const Array4<Real> mm(slopes, ncomp*AMREX_SPACEDIM);  // min and max

    for             (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cmn = u(i,j,k,nu);
                    Real cmx = cmn;
                    for         (int koff = -1; koff <= 1; ++koff) {
                        for     (int joff = -1; joff <= 1; ++joff) {
                            for (int ioff = -1; ioff <= 1; ++ioff) {
                                cmn = amrex::min(cmn,u(i+ioff,j+joff,k+koff,nu));
                                cmx = amrex::max(cmx,u(i+ioff,j+joff,k+koff,nu));
                            }
                        }
                    }
                    mm(i,j,k,n      ) = cmn - u(i,j,k,nu);
                    mm(i,j,k,n+ncomp) = cmx - u(i,j,k,nu);
                }
            }
        }

        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = Real(2.)*(u(i+1,j,k,nu)-u(i  ,j,k,nu));
                    Real back = Real(2.)*(u(i  ,j,k,nu)-u(i-1,j,k,nu));
                    Real slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = Real(2.)*(u(i,j+1,k,nu)-u(i,j  ,k,nu));
                    back = Real(2.)*(u(i,j  ,k,nu)-u(i,j-1,k,nu));
                    slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n+ncomp) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = Real(2.)*(u(i,j,k+1,nu)-u(i,j,k  ,nu));
                    back = Real(2.)*(u(i,j,k  ,nu)-u(i,j,k-1,nu));
                    slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                    slopes(i,j,k,n+ncomp*2) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cellconslin_fine_alpha (Box const& bx, Array4<Real> const& alpha,
                        Array4<Real const> const& slopes, const int ncomp,
                        Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Array4<Real const> mm(slopes, ncomp*AMREX_SPACEDIM);  // min and max

    const auto vlo  = amrex::lbound(alpha);
    const auto vlen = amrex::length(alpha);
    Real const* AMREX_RESTRICT xoff = voff;
    Real const* AMREX_RESTRICT yoff = voff + vlen.x;
    Real const* AMREX_RESTRICT zoff = voff + (vlen.x+vlen.y);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    const Real dummy_fine = xoff[i-vlo.x]*slopes(ic,jc,kc,n)
                        +                   yoff[j-vlo.y]*slopes(ic,jc,kc,n+ncomp)
                        +                   zoff[k-vlo.z]*slopes(ic,jc,kc,n+ncomp*2);

                    if (dummy_fine > mm(ic,jc,kc,n+ncomp) && dummy_fine != Real(0.)) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n+ncomp) / dummy_fine;
                    } else if (dummy_fine < mm(ic,jc,kc,n) && dummy_fine != Real(0.)) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n) / dummy_fine;
                    } else {
                        alpha(i,j,k,n) = Real(1.);
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cellconslin_slopes_mmlim (Box const& bx, Array4<Real> const& slopes,
                          Array4<Real const> const& alpha,
                          const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kk = k*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jj = j*ratio[1];
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ii = i*ratio[0];
                    Real a = Real(1.);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pcinterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp);
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

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/ratio[0];
    const Real ry = Real(1.)/ratio[1];
    const Real rz = Real(1.)/ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    T dx00 = u(i+1,j,k,nu) - u(i,j,k,nu);
                    T d0x0 = u(i,j+1,k,nu) - u(i,j,k,nu);
                    T d00x = u(i,j,k+1,nu) - u(i,j,k,nu);

                    T dx10 = u(i+1,j+1,k,nu) - u(i,j+1,k,nu);
                    T dx01 = u(i+1,j,k+1,nu) - u(i,j,k+1,nu);
                    T d0x1 = u(i,j+1,k+1,nu) - u(i,j,k+1,nu);

                    T dx11 = u(i+1,j+1,k+1,nu) - u(i,j+1,k+1,nu);

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

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_interp (Box const& bx, Array4<T> const& fine, const int fcomp, const int ncomp,
                  Array4<T const> const& slope, Array4<T const> const& crse,
                  const int ccomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto chi = amrex::ubound(slope);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::min(amrex::coarsen(k,ratio[2]),chi.z);
            const Real fz = k - kc*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
                const Real fy = j - jc*ratio[1];
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                    const Real fx = i - ic*ratio[0];
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp)
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

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_avgcrse_x (Box const& bx,
                   Array4<T> crse, Array4<T const> fine_val,
                   const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
/*  2D Copy (+ z loop added)

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int k = lo.z; k <= hi.z; ++k) {
            const int fk = k*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int fj = j*ratio[1];

                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int fi = i*ratio[0];

                    crse(i,j,0,nu) = ( fine_val(fi,fj,0,n) + fine_val(fi,fj+1,0,n) ) / 2.0_rt;
                }
            }
        }
    }
*/
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_avgcrse_y (Box const& bx,
                   Array4<T> crse, Array4<T const> fine_val,
                   const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_avgcrse_z (Box const& bx,
                   Array4<T> crse, Array4<T const> fine_val,
                   const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_ext (Box const& bx,
             GpuArray<Array4<T const>, AMREX_SPACEDIM> crse,
             GpuArray<Array4<T>, AMREX_SPACEDIM> fine,
             const int icomp, const int ncomp, IntVect const& ratio,
             const GpuArray<Real, AMREX_SPACEDIM> cellSize) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if ((cellSize[0] != cellSize[1]) || (cellSize[1] != cellSize[2]))
    {
        amrex::Abort( "FaceDivFree 3D stencil only built for dx=dy=dz." );
    }

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int k = lo.z; k <= hi.z; ++k) {
            const int fk = k*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int fj = j*ratio[1];

                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int fi = i*ratio[0];

                    // References to coarse values.
                    // Note the "origin" for each velocity
                    //     is the "negative" direction face of the current cell.
                    const Real& uW   = crse[0](i,   j,   k,   nu);
                    const Real& uWF  = crse[0](i,   j,   k+1, nu);
                    const Real& uWB  = crse[0](i,   j,   k-1, nu);
                    const Real& uNW  = crse[0](i,   j+1, k,   nu);
                    const Real& uNWF = crse[0](i,   j+1, k+1, nu);
                    const Real& uNWB = crse[0](i,   j+1, k-1, nu);
                    const Real& uSW  = crse[0](i,   j-1, k,   nu);
                    const Real& uSWF = crse[0](i,   j-1, k+1, nu);
                    const Real& uSWB = crse[0](i,   j-1, k-1, nu);

                    const Real& uE   = crse[0](i+1, j,   k,   nu);
                    const Real& uEF  = crse[0](i+1, j,   k+1, nu);
                    const Real& uEB  = crse[0](i+1, j,   k-1, nu);
                    const Real& uNE  = crse[0](i+1, j+1, k,   nu);
                    const Real& uNEF = crse[0](i+1, j+1, k+1, nu);
                    const Real& uNEB = crse[0](i+1, j+1, k-1, nu);
                    const Real& uSE  = crse[0](i+1, j-1, k,   nu);
                    const Real& uSEF = crse[0](i+1, j-1, k+1, nu);
                    const Real& uSEB = crse[0](i+1, j-1, k-1, nu);

                    const Real& vS   = crse[1](i,   j,   k,   nu);
                    const Real& vSF  = crse[1](i,   j,   k+1, nu);
                    const Real& vSB  = crse[1](i,   j,   k-1, nu);
                    const Real& vSW  = crse[1](i-1, j,   k,   nu);
                    const Real& vSWF = crse[1](i-1, j,   k+1, nu);
                    const Real& vSWB = crse[1](i-1, j,   k-1, nu);
                    const Real& vSE  = crse[1](i+1, j,   k,   nu);
                    const Real& vSEF = crse[1](i+1, j,   k+1, nu);
                    const Real& vSEB = crse[1](i+1, j,   k-1, nu);

                    const Real& vN   = crse[1](i,   j+1, k,   nu);
                    const Real& vNF  = crse[1](i,   j+1, k+1, nu);
                    const Real& vNB  = crse[1](i,   j+1, k-1, nu);
                    const Real& vNW  = crse[1](i-1, j+1, k,   nu);
                    const Real& vNWF = crse[1](i-1, j+1, k+1, nu);
                    const Real& vNWB = crse[1](i-1, j+1, k-1, nu);
                    const Real& vNE  = crse[1](i+1, j+1, k,   nu);
                    const Real& vNEF = crse[1](i+1, j+1, k+1, nu);
                    const Real& vNEB = crse[1](i+1, j+1, k-1, nu);

                    const Real& wB   = crse[2](i,     j,   k, nu);
                    const Real& wWB  = crse[2](i-1,   j,   k, nu);
                    const Real& wEB  = crse[2](i+1,   j,   k, nu);
                    const Real& wNB  = crse[2](i,   j+1,   k, nu);
                    const Real& wNWB = crse[2](i-1, j+1,   k, nu);
                    const Real& wNEB = crse[2](i+1, j+1,   k, nu);
                    const Real& wSB  = crse[2](i,   j-1,   k, nu);
                    const Real& wSWB = crse[2](i-1, j-1,   k, nu);
                    const Real& wSEB = crse[2](i+1, j-1,   k, nu);

                    const Real& wF   = crse[2](i,     j, k+1, nu);
                    const Real& wWF  = crse[2](i-1,   j, k+1, nu);
                    const Real& wEF  = crse[2](i+1,   j, k+1, nu);
                    const Real& wNF  = crse[2](i,   j+1, k+1, nu);
                    const Real& wNWF = crse[2](i-1, j+1, k+1, nu);
                    const Real& wNEF = crse[2](i+1, j+1, k+1, nu);
                    const Real& wSF  = crse[2](i,   j-1, k+1, nu);
                    const Real& wSWF = crse[2](i-1, j-1, k+1, nu);
                    const Real& wSEF = crse[2](i+1, j-1, k+1, nu);

                    // Calc fine faces that fall on coarse faces (exterior faces)
                    // All cells calc S and W and B exterior faces.
                    // N and E and F exterior faces only calced by edge cells.
                    fine[0](fi,  fj,  fk,  n) = (64.0_rt*uW + 8.0_rt*(-uNW + uSW - uWF + uWB) + uNWF - uNWB - uSWF + uSWB)/64.0_rt;
                    fine[0](fi,  fj+1,fk,  n) = (64.0_rt*uW + 8.0_rt*( uNW - uSW - uWF + uWB) - uNWF + uNWB + uSWF - uSWB)/64.0_rt;
                    fine[0](fi,  fj,  fk+1,n) = (64.0_rt*uW + 8.0_rt*(-uNW + uSW + uWF - uWB) - uNWF + uNWB + uSWF - uSWB)/64.0_rt;
                    fine[0](fi,  fj+1,fk+1,n) = (64.0_rt*uW + 8.0_rt*( uNW - uSW + uWF - uWB) + uNWF - uNWB - uSWF + uSWB)/64.0_rt;

                    fine[1](fi,  fj,  fk  ,n) = (64.0_rt*vS + 8.0_rt*( vSW - vSE - vSF + vSB) + vSWB - vSEB - vSWF + vSEF)/64.0_rt;
                    fine[1](fi+1,fj,  fk  ,n) = (64.0_rt*vS + 8.0_rt*(-vSW + vSE - vSF + vSB) - vSWB + vSEB + vSWF - vSEF)/64.0_rt;
                    fine[1](fi,  fj,  fk+1,n) = (64.0_rt*vS + 8.0_rt*( vSW - vSE + vSF - vSB) - vSWB + vSEB + vSWF - vSEF)/64.0_rt;
                    fine[1](fi+1,fj,  fk+1,n) = (64.0_rt*vS + 8.0_rt*(-vSW + vSE + vSF - vSB) + vSWB - vSEB - vSWF + vSEF)/64.0_rt;

                    fine[2](fi,  fj,  fk  ,n) = (64.0_rt*wB + 8.0_rt*(-wNB + wSB + wWB - wEB) - wNWB + wNEB + wSWB - wSEB)/64.0_rt;
                    fine[2](fi,  fj+1,fk  ,n) = (64.0_rt*wB + 8.0_rt*( wNB - wSB + wWB - wEB) + wNWB - wNEB - wSWB + wSEB)/64.0_rt;
                    fine[2](fi+1,fj,  fk  ,n) = (64.0_rt*wB + 8.0_rt*(-wNB + wSB - wWB + wEB) + wNWB - wNEB - wSWB + wSEB)/64.0_rt;
                    fine[2](fi+1,fj+1,fk  ,n) = (64.0_rt*wB + 8.0_rt*( wNB - wSB - wWB + wEB) - wNWB + wNEB + wSWB - wSEB)/64.0_rt;

                    if (i == hi.x)
                    {
                        fine[0](fi+2,fj,  fk,  n) = (64.0_rt*uE + 8.0_rt*(-uNE + uSE - uEF + uEB) + uNEF - uNEB - uSEF + uSEB)/64.0_rt;
                        fine[0](fi+2,fj+1,fk,  n) = (64.0_rt*uE + 8.0_rt*( uNE - uSE - uEF + uEB) - uNEF + uNEB + uSEF - uSEB)/64.0_rt;
                        fine[0](fi+2,fj,  fk+1,n) = (64.0_rt*uE + 8.0_rt*(-uNE + uSE + uEF - uEB) - uNEF + uNEB + uSEF - uSEB)/64.0_rt;
                        fine[0](fi+2,fj+1,fk+1,n) = (64.0_rt*uE + 8.0_rt*( uNE - uSE + uEF - uEB) + uNEF - uNEB - uSEF + uSEB)/64.0_rt;
                    }
                    if (j == hi.y)
                    {
                        fine[1](fi,  fj+2,fk  ,n) = (64.0_rt*vN + 8.0_rt*( vNW - vNE - vNF + vNB) + vNWB - vNEB - vNWF + vNEF)/64.0_rt;
                        fine[1](fi+1,fj+2,fk  ,n) = (64.0_rt*vN + 8.0_rt*(-vNW + vNE - vNF + vNB) - vNWB + vNEB + vNWF - vNEF)/64.0_rt;
                        fine[1](fi,  fj+2,fk+1,n) = (64.0_rt*vN + 8.0_rt*( vNW - vNE + vNF - vNB) - vNWB + vNEB + vNWF - vNEF)/64.0_rt;
                        fine[1](fi+1,fj+2,fk+1,n) = (64.0_rt*vN + 8.0_rt*(-vNW + vNE + vNF - vNB) + vNWB - vNEB - vNWF + vNEF)/64.0_rt;
                    }
                    if (k == hi.z)
                    {
                        fine[2](fi,  fj,  fk+2,n) = (64.0_rt*wF + 8.0_rt*(-wNF + wSF + wWF - wEF) - wNWF + wNEF + wSWF - wSEF)/64.0_rt;
                        fine[2](fi,  fj+1,fk+2,n) = (64.0_rt*wF + 8.0_rt*( wNF - wSF + wWF - wEF) + wNWF - wNEF - wSWF + wSEF)/64.0_rt;
                        fine[2](fi+1,fj,  fk+2,n) = (64.0_rt*wF + 8.0_rt*(-wNF + wSF - wWF + wEF) + wNWF - wNEF - wSWF + wSEF)/64.0_rt;
                        fine[2](fi+1,fj+1,fk+2,n) = (64.0_rt*wF + 8.0_rt*( wNF - wSF - wWF + wEF) - wNWF + wNEF + wSWF - wSEF)/64.0_rt;
                    }
                }
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (Box const& bx,
             GpuArray<Array4<T const>, AMREX_SPACEDIM> crse,
             GpuArray<Array4<T>, AMREX_SPACEDIM> fine,
             const int icomp, const int ncomp, IntVect const& ratio,
             const GpuArray<Real, AMREX_SPACEDIM> cellSize) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if ((cellSize[0] != cellSize[1]) || (cellSize[1] != cellSize[2]))
    {
        amrex::Abort( "FaceDivFree 3D stencil only built for dx=dy=dz." );
    }

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int k = lo.z; k <= hi.z; ++k) {
            const int fk = k*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int fj = j*ratio[1];

                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int fi = i*ratio[0];

                    // References to coarse values.
                    // Note the "origin" for each velocity
                    //     is the "negative" direction face of the current cell.
                    const Real& uW   = crse[0](i,   j,   k,   nu);
                    const Real& uWF  = crse[0](i,   j,   k+1, nu);
                    const Real& uWB  = crse[0](i,   j,   k-1, nu);
                    const Real& uNW  = crse[0](i,   j+1, k,   nu);
                    const Real& uNWF = crse[0](i,   j+1, k+1, nu);
                    const Real& uNWB = crse[0](i,   j+1, k-1, nu);
                    const Real& uSW  = crse[0](i,   j-1, k,   nu);
                    const Real& uSWF = crse[0](i,   j-1, k+1, nu);
                    const Real& uSWB = crse[0](i,   j-1, k-1, nu);

                    const Real& uE   = crse[0](i+1, j,   k,   nu);
                    const Real& uEF  = crse[0](i+1, j,   k+1, nu);
                    const Real& uEB  = crse[0](i+1, j,   k-1, nu);
                    const Real& uNE  = crse[0](i+1, j+1, k,   nu);
                    const Real& uNEF = crse[0](i+1, j+1, k+1, nu);
                    const Real& uNEB = crse[0](i+1, j+1, k-1, nu);
                    const Real& uSE  = crse[0](i+1, j-1, k,   nu);
                    const Real& uSEF = crse[0](i+1, j-1, k+1, nu);
                    const Real& uSEB = crse[0](i+1, j-1, k-1, nu);

                    const Real& vS   = crse[1](i,   j,   k,   nu);
                    const Real& vSF  = crse[1](i,   j,   k+1, nu);
                    const Real& vSB  = crse[1](i,   j,   k-1, nu);
                    const Real& vSW  = crse[1](i-1, j,   k,   nu);
                    const Real& vSWF = crse[1](i-1, j,   k+1, nu);
                    const Real& vSWB = crse[1](i-1, j,   k-1, nu);
                    const Real& vSE  = crse[1](i+1, j,   k,   nu);
                    const Real& vSEF = crse[1](i+1, j,   k+1, nu);
                    const Real& vSEB = crse[1](i+1, j,   k-1, nu);

                    const Real& vN   = crse[1](i,   j+1, k,   nu);
                    const Real& vNF  = crse[1](i,   j+1, k+1, nu);
                    const Real& vNB  = crse[1](i,   j+1, k-1, nu);
                    const Real& vNW  = crse[1](i-1, j+1, k,   nu);
                    const Real& vNWF = crse[1](i-1, j+1, k+1, nu);
                    const Real& vNWB = crse[1](i-1, j+1, k-1, nu);
                    const Real& vNE  = crse[1](i+1, j+1, k,   nu);
                    const Real& vNEF = crse[1](i+1, j+1, k+1, nu);
                    const Real& vNEB = crse[1](i+1, j+1, k-1, nu);

                    const Real& wB   = crse[2](i,     j,   k, nu);
                    const Real& wWB  = crse[2](i-1,   j,   k, nu);
                    const Real& wEB  = crse[2](i+1,   j,   k, nu);
                    const Real& wNB  = crse[2](i,   j+1,   k, nu);
                    const Real& wNWB = crse[2](i-1, j+1,   k, nu);
                    const Real& wNEB = crse[2](i+1, j+1,   k, nu);
                    const Real& wSB  = crse[2](i,   j-1,   k, nu);
                    const Real& wSWB = crse[2](i-1, j-1,   k, nu);
                    const Real& wSEB = crse[2](i+1, j-1,   k, nu);

                    const Real& wF   = crse[2](i,     j, k+1, nu);
                    const Real& wWF  = crse[2](i-1,   j, k+1, nu);
                    const Real& wEF  = crse[2](i+1,   j, k+1, nu);
                    const Real& wNF  = crse[2](i,   j+1, k+1, nu);
                    const Real& wNWF = crse[2](i-1, j+1, k+1, nu);
                    const Real& wNEF = crse[2](i+1, j+1, k+1, nu);
                    const Real& wSF  = crse[2](i,   j-1, k+1, nu);
                    const Real& wSWF = crse[2](i-1, j-1, k+1, nu);
                    const Real& wSEF = crse[2](i+1, j-1, k+1, nu);

                    // Calc interior fine faces from exterior values.
//                    const Real dxdy = cellSize[0]/cellSize[1];

                    fine[0](fi+1,fj  ,fk  ,n) = (1/256.0_rt)*(128.0_rt*(uW+uE)
                                                             + 16.0_rt*(  uWB +  uEB +  uSE +  uSW -  uEF -  uNE -  uNW -  uWF
                                                                        + vNE +  vSW -  vNW -  vSE +  wEF +  wWB -  wEB -  wWF)
                                                             +  2.0_rt*( uNEF + uNWF + uSEB + uSWB - uNEB - uNWB - uSEF - uSWF)
                                                             +         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF)
                                                             +         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[0](fi+1,fj+1,fk  ,n) = (1/256.0_rt)*(128.0_rt*(uW+uE)
                                                             + 16.0_rt*(  uWB +  uEB -  uSE -  uSW -  uEF +  uNE +  uNW -  uWF
                                                                        + vNE +  vSW -  vNW -  vSE +  wEF +  wWB -  wEB -  wWF)
                                                             -  2.0_rt*( uNEF + uNWF + uSEB + uSWB - uNEB - uNWB - uSEF - uSWF)
                                                             +         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF)
                                                             -         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[0](fi+1,fj  ,fk+1,n) = (1/256.0_rt)*(128.0_rt*(uW+uE)
                                                             + 16.0_rt*( -uWB -  uEB +  uSE +  uSW +  uEF -  uNE -  uNW +  uWF
                                                                        + vNE +  vSW -  vNW -  vSE +  wEF +  wWB -  wEB -  wWF)
                                                             -  2.0_rt*( uNEF + uNWF + uSEB + uSWB - uNEB - uNWB - uSEF - uSWF)
                                                             -         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF)
                                                             +         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[0](fi+1,fj+1,fk+1,n) = (1/256.0_rt)*(128.0_rt*(uW+uE)
                                                             + 16.0_rt*( -uWB -  uEB -  uSE -  uSW +  uEF +  uNE +  uNW +  uWF
                                                                        + vNE +  vSW -  vNW -  vSE +  wEF +  wWB -  wEB -  wWF)
                                                             +  2.0_rt*( uNEF + uNWF + uSEB + uSWB - uNEB - uNWB - uSEF - uSWF)
                                                             +         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF)
                                                             +         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[1](fi  ,fj+1,fk  ,n) = (1/256.0_rt)*(128.0_rt*(vN+vS)
                                                             + 16.0_rt*(  vNB +  vNW +  vSB +  vSW -  vNF -  vNE -  vSF -  vSE
                                                                        + uNE +  uSW -  uNW -  uSE +  wNF +  wSB -  wNB -  wSF)
                                                             +  2.0_rt*( vNEF + vNWB + vSEF + vSWB - vNEB - vNWF - vSEB - vSWF)
                                                             +         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             +         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[1](fi  ,fj+1,fk+1,n) = (1/256.0_rt)*(128.0_rt*(vN+vS)
                                                             + 16.0_rt*( -vNB +  vNW -  vSB +  vSW +  vNF -  vNE +  vSF -  vSE
                                                                        + uNE +  uSW -  uNW -  uSE +  wNF +  wSB -  wNB -  wSF)
                                                             -  2.0_rt*( vNEF + vNWB + vSEF + vSWB - vNEB - vNWF - vSEB - vSWF)
                                                             +         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             -         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[1](fi+1,fj+1,fk  ,n) = (1/256.0_rt)*(128.0_rt*(vN+vS)
                                                             + 16.0_rt*(  vNB -  vNW +  vSB -  vSW -  vNF +  vNE -  vSF +  vSE
                                                                        + uNE +  uSW -  uNW -  uSE +  wNF +  wSB -  wNB -  wSF)
                                                             -  2.0_rt*( vNEF + vNWB + vSEF + vSWB - vNEB - vNWF - vSEB - vSWF)
                                                             -         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             +         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[1](fi+1,fj+1,fk+1,n) = (1/256.0_rt)*(128.0_rt*(vN+vS)
                                                             + 16.0_rt*( -vNB -  vNW -  vSB -  vSW +  vNF +  vNE +  vSF +  vSE
                                                                        + uNE +  uSW -  uNW -  uSE +  wNF +  wSB -  wNB -  wSF)
                                                             +  2.0_rt*( vNEF + vNWB + vSEF + vSWB - vNEB - vNWF - vSEB - vSWF)
                                                             -         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             -         ( wNEB + wNWF + wSEF + wSWB - wNEF - wNWB - wSEB - wSWF));

                    fine[2](fi  ,fj  ,fk+1,n) = (1/256.0_rt)*(128.0_rt*(wF+wB)
                                                             + 16.0_rt*(  wSF +  wSB +  wWF +  wWB -  wEF -  wEB -  wNF -  wNB
                                                                        + uEF +  uWB -  uEB -  uWF +  vNF +  vSB -  vNB -  vSF)
                                                             +  2.0_rt*( wNEF + wNEB + wSWF + wSWB - wNWF - wNWB - wSEF - wSEB)
                                                             +         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             +         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF));

                    fine[2](fi  ,fj+1,fk+1,n) = (1/256.0_rt)*(128.0_rt*(wF+wB)
                                                             + 16.0_rt*( -wSF -  wSB +  wWF +  wWB -  wEF -  wEB +  wNF +  wNB
                                                                        + uEF +  uWB -  uEB -  uWF +  vNF +  vSB -  vNB -  vSF)
                                                             -  2.0_rt*( wNEF + wNEB + wSWF + wSWB - wNWF - wNWB - wSEF - wSEB)
                                                             +         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             -         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF));

                    fine[2](fi+1,fj  ,fk+1,n) = (1/256.0_rt)*(128.0_rt*(wF+wB)
                                                             + 16.0_rt*(  wSF +  wSB -  wWF -  wWB +  wEF +  wEB -  wNF -  wNB
                                                                        + uEF +  uWB -  uEB -  uWF +  vNF +  vSB -  vNB -  vSF)
                                                             -  2.0_rt*( wNEF + wNEB + wSWF + wSWB - wNWF - wNWB - wSEF - wSEB)
                                                             -         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             +         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF));

                    fine[2](fi+1,fj+1,fk+1,n) = (1/256.0_rt)*(128.0_rt*(wF+wB)
                                                             + 16.0_rt*( -wSF -  wSB -  wWF -  wWB +  wEF +  wEB +  wNF +  wNB
                                                                        + uEF +  uWB -  uEB -  uWF +  vNF +  vSB -  vNB -  vSF)
                                                             +  2.0_rt*( wNEF + wNEB + wSWF + wSWB - wNWF - wNWB - wSEF - wSEB)
                                                             -         ( uNEB + uNWF + uSEF + uSWB - uNEF - uNWB - uSEB - uSWF)
                                                             -         ( vNEB + vNWF + vSEF + vSWB - vNEF - vNWB - vSEB - vSWF));
                }
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    int ii = amrex::coarsen(i,ratio[0]);
    int jj = amrex::coarsen(j,ratio[1]);
    int kk = amrex::coarsen(k,ratio[2]);
    Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
    fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii+1,jj,kk,n);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    int ii = amrex::coarsen(i,ratio[0]);
    int jj = amrex::coarsen(j,ratio[1]);
    int kk = amrex::coarsen(k,ratio[2]);
    Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/ratio[1]);
    fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj+1,kk,n);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_z (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    int ii = amrex::coarsen(i,ratio[0]);
    int jj = amrex::coarsen(j,ratio[1]);
    int kk = amrex::coarsen(k,ratio[2]);
    Real const w = static_cast<Real>(k-kk*ratio[2]) * (Real(1.)/ratio[2]);
    fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj,kk+1,n);
}


}  // namespace amrex


#endif
