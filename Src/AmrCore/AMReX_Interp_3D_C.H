#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_
#include <AMReX_Config.H>

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
facediv_face_interp (int ci, int cj, int ck,
                     int nc, int nf, int idir,
                     Array4<T const> const& crse,
                     Array4<T> const& fine,
                     Array4<int const> const& mask,
                     IntVect const& ratio) noexcept
{

    if (mask) {
        if (!mask(ci, cj, ck, nc))
            { return; }
    }

    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];
    const int fk = ck*ratio[2];

    switch (idir) {
        case 0:
        {
            const Real ll = crse(ci, cj-1, ck-1, nc);
            const Real cl = crse(ci, cj-1, ck,   nc);
            const Real rl = crse(ci, cj-1, ck+1, nc);

            const Real lc = crse(ci, cj  , ck-1, nc);
            const Real cc = crse(ci, cj  , ck,   nc);
            const Real rc = crse(ci, cj  , ck+1, nc);

            const Real lu = crse(ci, cj+1, ck-1, nc);
            const Real cu = crse(ci, cj+1, ck,   nc);
            const Real ru = crse(ci, cj+1, ck+1, nc);

            fine(fi,fj  ,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi,fj  ,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi,fj+1,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi,fj+1,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
        case 1:
        {
            const Real ll = crse(ci-1, cj, ck-1, nc);
            const Real cl = crse(ci  , cj, ck-1, nc);
            const Real rl = crse(ci+1, cj, ck-1, nc);

            const Real lc = crse(ci-1, cj, ck  , nc);
            const Real cc = crse(ci  , cj, ck  , nc);
            const Real rc = crse(ci+1, cj, ck  , nc);

            const Real lu = crse(ci-1, cj, ck+1, nc);
            const Real cu = crse(ci  , cj, ck+1, nc);
            const Real ru = crse(ci+1, cj, ck+1, nc);

            fine(fi  ,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi+1,fj,fk  ,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi  ,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi+1,fj,fk+1,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
        case 2:
        {
            const Real ll = crse(ci-1, cj-1, ck, nc);
            const Real cl = crse(ci  , cj-1, ck, nc);
            const Real rl = crse(ci+1, cj-1, ck, nc);

            const Real lc = crse(ci-1, cj  , ck, nc);
            const Real cc = crse(ci  , cj  , ck, nc);
            const Real rc = crse(ci+1, cj  , ck, nc);

            const Real lu = crse(ci-1, cj+1, ck, nc);
            const Real cu = crse(ci  , cj+1, ck, nc);
            const Real ru = crse(ci+1, cj+1, ck, nc);

            fine(fi  ,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+lc-cu-rc) + (ll+ru-lu-rl) );
            fine(fi+1,fj  ,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cl+rc-cu-lc) + (lu+rl-ll-ru) );
            fine(fi  ,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+lc-cl-rc) + (lu+rl-ll-ru) );
            fine(fi+1,fj+1,fk,nf) = Real(1.)/Real(64.)*( 64*cc + 8*(cu+rc-cl-lc) + (ll+ru-lu-rl) );

            break;
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int ci, int cj, int ck, int nf,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& fine,
             IntVect const& ratio,
             GpuArray<Real, AMREX_SPACEDIM> const& cellSize) noexcept
{
    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];
    const int fk = ck*ratio[1];

    // References to fine exterior values needed for interior calculation.
    const Real u000 = fine[0](fi,   fj  , fk  , nf);
    const Real u200 = fine[0](fi+2, fj  , fk  , nf);
    const Real u010 = fine[0](fi,   fj+1, fk  , nf);
    const Real u210 = fine[0](fi+2, fj+1, fk  , nf);
    const Real u001 = fine[0](fi,   fj  , fk+1, nf);
    const Real u201 = fine[0](fi+2, fj  , fk+1, nf);
    const Real u011 = fine[0](fi,   fj+1, fk+1, nf);
    const Real u211 = fine[0](fi+2, fj+1, fk+1, nf);

    const Real v000 = fine[1](fi  , fj  , fk  , nf);
    const Real v020 = fine[1](fi  , fj+2, fk  , nf);
    const Real v100 = fine[1](fi+1, fj  , fk  , nf);
    const Real v120 = fine[1](fi+1, fj+2, fk  , nf);
    const Real v001 = fine[1](fi  , fj  , fk+1, nf);
    const Real v021 = fine[1](fi  , fj+2, fk+1, nf);
    const Real v101 = fine[1](fi+1, fj  , fk+1, nf);
    const Real v121 = fine[1](fi+1, fj+2, fk+1, nf);

    const Real w000 = fine[2](fi  , fj  , fk  , nf);
    const Real w002 = fine[2](fi  , fj  , fk+2, nf);
    const Real w100 = fine[2](fi+1, fj  , fk  , nf);
    const Real w102 = fine[2](fi+1, fj  , fk+2, nf);
    const Real w010 = fine[2](fi  , fj+1, fk  , nf);
    const Real w012 = fine[2](fi  , fj+1, fk+2, nf);
    const Real w110 = fine[2](fi+1, fj+1, fk  , nf);
    const Real w112 = fine[2](fi+1, fj+1, fk+2, nf);

    const Real dx = cellSize[0];
    const Real dy = cellSize[1];
    const Real dz = cellSize[2];

    const Real dx3 = dx*dx*dx;
    const Real dy3 = dy*dy*dy;
    const Real dz3 = dz*dz*dz;

    // aspbs = a squared plus b squared
    const Real xspys = dx*dx + dy*dy;
    const Real yspzs = dy*dy + dz*dz;
    const Real zspxs = dz*dz + dx*dx;

    fine[0](fi+1, fj  , fk  , nf) = Real(0.5)*(u000+u200)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  +                dx3/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w000+w102-w002-w100)
                                  +                dx3/(8*dz*xspys)*(w010+w112-w012-w110);

    fine[0](fi+1, fj+1, fk  , nf) = Real(0.5)*(u010+u210)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  +                dx3/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w010+w112-w012-w110)
                                  +                dx3/(8*dz*xspys)*(w000+w102-w100-w002);

    fine[0](fi+1, fj  , fk+1, nf) = Real(0.5)*(u001+u201)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  +                dx3/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w000+w102-w002-w100)
                                  +                dx3/(8*dz*xspys)*(w010+w112-w012-w110);

    fine[0](fi+1, fj+1, fk+1, nf) = Real(0.5)*(u011+u211)
                                  + dx*(2*dz*dz+dx*dx)/(8*dy*zspxs)*(v001+v121-v021-v101)
                                  +                dx3/(8*dy*zspxs)*(v000+v120-v020-v100)
                                  + dx*(2*dy*dy+dx*dx)/(8*dz*xspys)*(w010+w112-w012-w110)
                                  +                dx3/(8*dz*xspys)*(w000+w102-w002-w100);

    fine[1](fi  , fj+1, fk  , nf) = Real(0.5)*(v000+v020)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  +                dy3/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w000+w012-w002-w010)
                                  +                dy3/(8*dz*xspys)*(w100+w112-w102-w110);

    fine[1](fi+1, fj+1, fk  , nf) = Real(0.5)*(v001+v021)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  +                dy3/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w100+w112-w102-w110)
                                  +                dy3/(8*dz*xspys)*(w000+w012-w002-w010);

    fine[1](fi  , fj+1, fk+1, nf) = Real(0.5)*(v100+v120)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  +                dy3/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w000+w012-w002-w010)
                                  +                dy3/(8*dz*xspys)*(w100+w112-w102-w110);

    fine[1](fi+1, fj+1, fk+1, nf) = Real(0.5)*(v101+v121)
                                  + dy*(2*dz*dz+dy*dy)/(8*dx*yspzs)*(u001+u211-u011-u201)
                                  +                dy3/(8*dx*yspzs)*(u000+u210-u010-u200)
                                  + dy*(2*dx*dx+dy*dy)/(8*dz*xspys)*(w100+w112-w102-w110)
                                  +                dy3/(8*dz*xspys)*(w000+w012-w002-w010);

    fine[2](fi  , fj  , fk+1, nf) = Real(0.5)*(w000+w002)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  +                dz3/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v000+v021-v001-v020)
                                  +                dz3/(8*dy*zspxs)*(v100+v121-v101-v120);

    fine[2](fi  , fj+1, fk+1, nf) = Real(0.5)*(w010+w012)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  +                dz3/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);

    fine[2](fi+1, fj  , fk+1, nf) = Real(0.5)*(w100+w102)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v000+v021-v001-v020)
                                  +                dz3/(8*dy*zspxs)*(v100+v121-v101-v120);

    fine[2](fi+1, fj+1, fk+1, nf) = Real(0.5)*(w110+w112)
                                  + dz*(2*dy*dy+dz*dz)/(8*dx*yspzs)*(u010+u211-u011-u210)
                                  +                dz3/(8*dx*yspzs)*(u000+u201-u001-u200)
                                  + dz*(2*dx*dx+dz*dz)/(8*dy*zspxs)*(v100+v121-v101-v120)
                                  +                dz3/(8*dy*zspxs)*(v000+v021-v001-v020);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (i-ii*ratio[0] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii+1,jj,kk,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (j-jj*ratio[1] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/ratio[1]);
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj+1,kk,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_z (int i, int j, int k, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    const int kk = amrex::coarsen(k,ratio[2]);
    if (k-kk*ratio[2] == 0) {
        fine(i,j,k,n) = crse(ii,jj,kk,n);
    } else {
        Real const w = static_cast<Real>(k-kk*ratio[2]) * (Real(1.)/ratio[2]);
        fine(i,j,k,n) = (Real(1.)-w) * crse(ii,jj,kk,n) + w * crse(ii,jj,kk+1,n);
    }
}

}  // namespace amrex


#endif
