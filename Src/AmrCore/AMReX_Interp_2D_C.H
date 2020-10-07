#ifndef AMREX_INTERP_2D_C_H_
#define AMREX_INTERP_2D_C_H_

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
    const int nxy = flen.x + flen.y;
    Vector<Real> voff(nxy);

    Vector<Real> fvc, cvc;
    fgeom.GetEdgeVolCoord(fvc,fbx,0);
    cgeom.GetEdgeVolCoord(cvc,cbx,0);

    const auto& clo = amrex::lbound(cbx);

    Real* AMREX_RESTRICT xoff = &voff[0];
    AMREX_PRAGMA_SIMD
    for (int i = flo.x; i <= fhi.x; ++i) {
        const int ic = amrex::coarsen(i, ratio[0]);
        const int ii = i - flo.x;
        const int iic = ic - clo.x;
        const Real fcen = Real(0.5)*(fvc[ii ]+fvc[ii +1]);
        const Real ccen = Real(0.5)*(cvc[iic]+cvc[iic+1]);
        xoff[ii] = (fcen-ccen)/(cvc[iic+1]-cvc[iic]);
    }

    fgeom.GetEdgeVolCoord(fvc,fbx,1);
    cgeom.GetEdgeVolCoord(cvc,cbx,1);

    Real* AMREX_RESTRICT yoff = &voff[flen.x];
    AMREX_PRAGMA_SIMD
    for (int j = flo.y; j <= fhi.y; ++j) {
        const int jc = amrex::coarsen(j, ratio[1]);
        const int jj = j - flo.y;
        const int jjc = jc - clo.y;
        const Real fcen = Real(0.5)*(fvc[jj ]+fvc[jj +1]);
        const Real ccen = Real(0.5)*(cvc[jjc]+cvc[jjc+1]);
        yoff[jj] = (fcen-ccen)/(cvc[jjc+1]-cvc[jjc]);
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
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            slopes(i,j,0,ns      ) = Real(0.5)*(u(i+1,j,0,nu)-u(i-1,j,0,nu));
            slopes(i,j,0,ns+ncomp) = Real(0.5)*(u(i,j+1,0,nu)-u(i,j-1,0,nu));
        }
    }

    if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
    {
        const int i = slo.x;
        if (shi.x-slo.x >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                slopes(i,j,0,ns) = -Real(16./15.)*u(i-1,j,0,nu) + Real(0.5)*u(i,j,0,nu)
                    + Real(2./3.)*u(i+1,j,0,nu) - Real(0.1)*u(i+2,j,0,nu);
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                slopes(i,j,0,ns) = Real(0.25)*(u(i+1,j,0,nu)+Real(5.)*u(i,j,0,nu)-Real(6.)*u(i-1,j,0,nu));
            }
        }
    }

    if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
    {
        const int i = shi.x;
        if (shi.x-slo.x >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                slopes(i,j,0,ns) = Real(16./15.)*u(i+1,j,0,nu) - Real(0.5)*u(i,j,0,nu)
                    - Real(2./3.)*u(i-1,j,0,nu) + Real(0.1)*u(i-2,j,0,nu);
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                slopes(i,j,0,ns) = -Real(0.25)*(u(i-1,j,0,nu)+Real(5.)*u(i,j,0,nu)-Real(6.)*u(i+1,j,0,nu));
            }
        }
    }

    if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
    {
        const int j = slo.y;
        if (shi.y-slo.y >= 1) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,0,ns+ncomp) = -Real(16./15.)*u(i,j-1,0,nu) + Real(0.5)*u(i,j,0,nu)
                    + Real(2./3.)*u(i,j+1,0,nu) - Real(0.1)*u(i,j+2,0,nu);
            }
        } else {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,0,ns+ncomp) = Real(0.25)*(u(i,j+1,0,nu)+Real(5.)*u(i,j,0,nu)-Real(6.)*u(i,j-1,0,nu));
            }
        }
    }

    if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
    {
        const int j = shi.y;
        if (shi.y-slo.y >= 1) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,0,ns+ncomp) = Real(16./15.)*u(i,j+1,0,nu) - Real(0.5)*u(i,j,0,nu)
                    - Real(2./3.)*u(i,j-1,0,nu) + Real(0.1)*u(i,j-2,0,nu);
            }
        } else {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,0,ns+ncomp) = -Real(0.25)*(u(i,j-1,0,nu)+Real(5.)*u(i,j,0,nu)-Real(6.)*u(i,j+1,0,nu));
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

    const Array4<Real> sf(slopes,ncomp*AMREX_SPACEDIM); // slope factor

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            sf(i,j,0,0) = Real(1.);
            sf(i,j,0,1) = Real(1.);
        }
    }

    for (int n = 0; n < ncomp; ++n)
    {
        const int nu = n + icomp;
        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real cen  = slopes(i,j,0,n);
                Real forw = Real(2.)*(u(i+1,j,0,nu)-u(i  ,j,0,nu));
                Real back = Real(2.)*(u(i  ,j,0,nu)-u(i-1,j,0,nu));
                Real slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                slopes(i,j,0,n) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                if (cen != Real(0.)) {
                    sf(i,j,0,0) = amrex::min(sf(i,j,0,0), slopes(i,j,0,n)/cen);
                } else {
                    sf(i,j,0,0) = Real(0.);
                }

                cen  = slopes(i,j,0,n+ncomp);
                forw = Real(2.)*(u(i,j+1,0,nu)-u(i,j  ,0,nu));
                back = Real(2.)*(u(i,j  ,0,nu)-u(i,j-1,0,nu));
                slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                slopes(i,j,0,n+ncomp) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
                if (cen != Real(0.)) {
                    sf(i,j,0,1) = amrex::min(sf(i,j,0,1), slopes(i,j,0,n+ncomp)/cen);
                } else {
                    sf(i,j,0,1) = Real(0.);
                }
            }
        }
    }

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,0,n      ) *= sf(i,j,0,0);
                slopes(i,j,0,n+ncomp) *= sf(i,j,0,1);
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

    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp)
                    + xoff[i-vlo.x] * slopes(ic,jc,0,n)
                    + yoff[j-vlo.y] * slopes(ic,jc,0,n+ncomp);
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

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real cmn = u(i,j,0,nu);
                Real cmx = cmn;
                for     (int joff = -1; joff <= 1; ++joff) {
                    for (int ioff = -1; ioff <= 1; ++ioff) {
                        cmn = amrex::min(cmn,u(i+ioff,j+joff,0,nu));
                        cmx = amrex::max(cmx,u(i+ioff,j+joff,0,nu));
                    }
                }
                mm(i,j,0,n      ) = cmn - u(i,j,0,nu);
                mm(i,j,0,n+ncomp) = cmx - u(i,j,0,nu);
            }
        }

        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real cen  = slopes(i,j,0,n);
                Real forw = Real(2.)*(u(i+1,j,0,nu)-u(i  ,j,0,nu));
                Real back = Real(2.)*(u(i  ,j,0,nu)-u(i-1,j,0,nu));
                Real slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                slopes(i,j,0,n) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));

                cen  = slopes(i,j,0,n+ncomp);
                forw = Real(2.)*(u(i,j+1,0,nu)-u(i,j  ,0,nu));
                back = Real(2.)*(u(i,j  ,0,nu)-u(i,j-1,0,nu));
                slp = (forw*back >= Real(0.)) ? amrex::min(amrex::Math::abs(forw),amrex::Math::abs(back)) : Real(0.);
                slopes(i,j,0,n+ncomp) = amrex::Math::copysign(Real(1.),cen)*amrex::min(slp,amrex::Math::abs(cen));
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

    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                const Real dummy_fine = xoff[i-vlo.x]*slopes(ic,jc,0,n)
                    +                   yoff[j-vlo.y]*slopes(ic,jc,0,n+ncomp);

                if (dummy_fine > mm(ic,jc,0,n+ncomp) && dummy_fine != Real(0.)) {
                    alpha(i,j,0,n) = mm(ic,jc,0,n+ncomp) / dummy_fine;
                } else if (dummy_fine < mm(ic,jc,0,n) && dummy_fine != Real(0.)) {
                    alpha(i,j,0,n) = mm(ic,jc,0,n) / dummy_fine;
                } else {
                    alpha(i,j,0,n) = Real(1.);
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
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jj = j*ratio[1];
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ii = i*ratio[0];
                Real a = Real(1.);
                for     (int joff = 0; joff < ratio[1]; ++joff) {
                    for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                        a = amrex::min(a, alpha(ii+ioff,jj+joff,0,n));
                    }
                }
                slopes(i,j,0,n      ) *= a;
                slopes(i,j,0,n+ncomp) *= a;
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
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
            }
        }
    }
}

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int ixy  = 2;
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

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                T dx0 = u(i+1,j,0,n+icomp) - u(i,j,0,n+icomp);
                T d0x = u(i,j+1,0,n+icomp) - u(i,j,0,n+icomp);
                T dx1 = u(i+1,j+1,0,n+icomp) - u(i,j+1,0,n+icomp);

                slope(i,j,0,n+ncomp*ix ) = rx*dx0;
                slope(i,j,0,n+ncomp*iy ) = ry*d0x;
                slope(i,j,0,n+ncomp*ixy) = rx*ry*(dx1 - dx0);
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
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
            const Real fy = j - jc*ratio[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                const Real fx = i - ic*ratio[0];
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp)
                    + fx*slope(ic,jc,0,n+ncomp*ix)
                    + fy*slope(ic,jc,0,n+ncomp*iy)
                    + fx*fy*slope(ic,jc,0,n+ncomp*ixy);
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
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
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

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_avgcrse_y (Box const& bx,
                   Array4<T> crse, Array4<T const> fine_val,
                   const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int j = lo.y; j <= hi.y; ++j) {
            const int fj = j*ratio[1];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int fi = i*ratio[0];
                crse(i,j,0,nu) = ( fine_val(fi,fj,0,n) + fine_val(fi+1,fj,0,n) ) / 2.0_rt;
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_ext (Box const& bx,
             GpuArray<Array4<T const>, AMREX_SPACEDIM> crse,
             GpuArray<Array4<T>, AMREX_SPACEDIM> fine,
             const int icomp, const int ncomp, IntVect const& ratio,
             const GpuArray<Real, AMREX_SPACEDIM> /*cellSize*/) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int j = lo.y; j <= hi.y; ++j) {
            const int fj = j*ratio[1];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int fi = i*ratio[0];

                // References to coarse values.
                const Real& uW  = crse[0](i,   j,   0, nu);
                const Real& uNW = crse[0](i,   j+1, 0, nu);
                const Real& uSW = crse[0](i,   j-1, 0, nu);
                const Real& uE  = crse[0](i+1, j,   0, nu);
                const Real& uNE = crse[0](i+1, j+1, 0, nu);
                const Real& uSE = crse[0](i+1, j-1, 0, nu);

                const Real& vS  = crse[1](i,   j,   0, nu);
                const Real& vSW = crse[1](i-1, j,   0, nu);
                const Real& vSE = crse[1](i+1, j,   0, nu);
                const Real& vN  = crse[1](i,   j+1, 0, nu);
                const Real& vNW = crse[1](i-1, j+1, 0, nu);
                const Real& vNE = crse[1](i+1, j+1, 0, nu);

                // Calc fine faces on negative exterior faces of coarse cells (W and S).
                fine[0](fi,  fj  ,0,n) = (8.0_rt*uW - uNW + uSW)/8.0_rt;
                fine[0](fi,  fj+1,0,n) = (8.0_rt*uW + uNW - uSW)/8.0_rt;
                fine[1](fi,  fj  ,0,n) = (8.0_rt*vS + vSW - vSE)/8.0_rt;
                fine[1](fi+1,fj  ,0,n) = (8.0_rt*vS - vSW + vSE)/8.0_rt;

                // Calc fine faces on positive exterior faces of coarse cells (E and N).
                // Only done by cells touching those boundaries.
                if (i == hi.x)
                {
                    fine[0](fi+2,fj  ,0,n) = (8.0_rt*uE - uNE + uSE)/8.0_rt;
                    fine[0](fi+2,fj+1,0,n) = (8.0_rt*uE + uNE - uSE)/8.0_rt;
                }
                if (j == hi.y)
                {
                    fine[1](fi,  fj+2,0,n) = (8.0_rt*vN + vNW - vNE)/8.0_rt;
                    fine[1](fi+1,fj+2,0,n) = (8.0_rt*vN - vNW + vNE)/8.0_rt;
                }
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (Box const& bx,
             GpuArray<Array4<T const>, AMREX_SPACEDIM> /*crse*/,
             GpuArray<Array4<T>, AMREX_SPACEDIM> fine,
             const int icomp, const int ncomp, IntVect const& ratio,
             const GpuArray<Real, AMREX_SPACEDIM> cellSize) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for (int j = lo.y; j <= hi.y; ++j) {
            const int fj = j*ratio[1];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int fi = i*ratio[0];

                // References to fine exterior values.
                const Real& umm = fine[0](fi,   fj,   0, nu);
                const Real& ump = fine[0](fi,   fj+1, 0, nu);
                const Real& upm = fine[0](fi+2, fj,   0, nu);
                const Real& upp = fine[0](fi+2, fj+1, 0, nu);

                const Real& vmm = fine[1](fi,   fj,   0, nu);
                const Real& vmp = fine[1](fi+1, fj,   0, nu);
                const Real& vpm = fine[1](fi,   fj+2, 0, nu);
                const Real& vpp = fine[1](fi+1, fj+2, 0, nu);

                const Real dxdy = cellSize[0]/cellSize[1];
                const Real x_corr =   (dxdy/4.0_rt) * (vpp+vmm-vmp-vpm);
                const Real y_corr = 1/(dxdy*4.0_rt) * (upp+umm-ump-upm);

                // Calc fine faces on interior of coarse cells.
                fine[0](fi+1,fj  ,0,n) = (umm+upm)/2.0_rt + x_corr; 
                fine[0](fi+1,fj+1,0,n) = (ump+upp)/2.0_rt + x_corr; 
                fine[1](fi,  fj+1,0,n) = (vmm+vpm)/2.0_rt + y_corr;
                fine[1](fi+1,fj+1,0,n) = (vmp+vpp)/2.0_rt + y_corr;
            }
        }
    }
}



template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    int ii = amrex::coarsen(i,ratio[0]);
    int jj = amrex::coarsen(j,ratio[1]);
    Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
    fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii+1,jj,0,n);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    int ii = amrex::coarsen(i,ratio[0]);
    int jj = amrex::coarsen(j,ratio[1]);
    Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/ratio[1]);
    fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii,jj+1,0,n);
}

}

#endif
