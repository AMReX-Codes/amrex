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
        const Real fcen = 0.5*(fvc[ii ]+fvc[ii +1]);
        const Real ccen = 0.5*(cvc[iic]+cvc[iic+1]);
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
        const Real fcen = 0.5*(fvc[jj ]+fvc[jj +1]);
        const Real ccen = 0.5*(cvc[jjc]+cvc[jjc+1]);
        yoff[jj] = (fcen-ccen)/(cvc[jjc+1]-cvc[jjc]);
    }

    return voff;
}

namespace {

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
compute_slopes (const Dim3& len, const Dim3& lo, const Dim3& hi, FabView<Real> const& slopes,
                FabView<const Real> const& u, Box const& sbox, const Dim3& slo, const Dim3& shi,
                const Dim3& slen, BCRec const& bc, const int n, const int ncomp) noexcept
{
    for     (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            slopes(i,j,0,n      ) = 0.5*(u(i+1,j,0,n)-u(i-1,j,0,n));
            slopes(i,j,0,n+ncomp) = 0.5*(u(i,j+1,0,n)-u(i,j-1,0,n));
        }
    }

    if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
    {
        const int i = 0;
        if (slen.x >= 2) {
            for     (int j = 0; j < len.y; ++j) {
                slopes(i,j,0,n) = -(16./15.)*u(i-1,j,0,n) + 0.5*u(i,j,0,n)
                    + (2./3.)*u(i+1,j,0,n) - 0.1*u(i+2,j,0,n);
            }
        } else {
            for     (int j = 0; j < len.y; ++j) {
                slopes(i,j,0,n) = 0.25*(u(i+1,j,0,n)+5.*u(i,j,0,n)-6.*u(i-1,j,0,n));
            }
        }
    }

    if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
    {
        const int i = len.x-1;
        if (slen.x >= 2) {
            for     (int j = 0; j < len.y; ++j) {
                slopes(i,j,0,n) = (16./15.)*u(i+1,j,0,n) - 0.5*u(i,j,0,n)
                    - (2./3.)*u(i-1,j,0,n) + 0.1*u(i-2,j,0,n);
            }
        } else {
            for     (int j = 0; j < len.y; ++j) {
                slopes(i,j,0,n) = -0.25*(u(i-1,j,0,n)+5.*u(i,j,0,n)-6.*u(i+1,j,0,n));
            }
        }
    }

    if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
    {
        const int j = 0;
        if (slen.y >= 2) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,0,n+ncomp) = -(16./15.)*u(i,j-1,0,n) + 0.5*u(i,j,0,n)
                    + (2./3.)*u(i,j+1,0,n) - 0.1*u(i,j+2,0,n);
            }
        } else {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,0,n+ncomp) = 0.25*(u(i,j+1,0,n)+5.*u(i,j,0,n)-6.*u(i,j-1,0,n));
            }
        }
    }

    if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
    {
        const int j = len.y-1;
        if (slen.y >= 2) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,0,n+ncomp) = (16./15.)*u(i,j+1,0,n) - 0.5*u(i,j,0,n)
                    - (2./3.)*u(i,j-1,0,n) + 0.1*u(i,j-2,0,n);
            }
        } else {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,0,n+ncomp) = -0.25*(u(i,j-1,0,n)+5.*u(i,j,0,n)-6.*u(i,j+1,0,n));
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

    for     (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            sf(i,j,0,0) = 1.0;
            sf(i,j,0,1) = 1.0;
        }
    }

    for (int n = 0; n < ncomp; ++n)
    {
        compute_slopes(len, lo, hi, slopes, u, sbox, slo, shi, slen, bcr[n], n, ncomp);

        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real cen  = slopes(i,j,0,n);
                Real forw = 2.0*(u(i+1,j,0,n)-u(i  ,j,0,n));
                Real back = 2.0*(u(i  ,j,0,n)-u(i-1,j,0,n));
                Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                slopes(i,j,0,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                if (cen != 0.0) {
                    sf(i,j,0,0) = amrex::min(sf(i,j,0,0), slopes(i,j,0,n)/cen);
                } else {
                    sf(i,j,0,0) = 0.0;
                }

                cen  = slopes(i,j,0,n+ncomp);
                forw = 2.0*(u(i,j+1,0,n)-u(i,j  ,0,n));
                back = 2.0*(u(i,j  ,0,n)-u(i,j-1,0,n));
                slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                slopes(i,j,0,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                if (cen != 0.0) {
                    sf(i,j,0,1) = amrex::min(sf(i,j,0,1), slopes(i,j,0,n+ncomp)/cen);
                } else {
                    sf(i,j,0,1) = 0.0;
                }
            }
        }
    }

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                slopes(i,j,0,n      ) *= sf(i,j,0,0);
                slopes(i,j,0,n+ncomp) *= sf(i,j,0,1);
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

    const auto clo = amrex::coarsen(lo,ratio);
    const auto slopes = slopesfab.view(clo);
    const auto crse = crsefab.view(clo,ccomp);

    Box vbox = slopesfab.box();
    vbox.refine(ratio);
    const auto vlo  = amrex::lbound(vbox);
    const auto vlen = amrex::length(vbox);
    Real const* AMREX_RESTRICT xoff = &voff[lo.x-vlo.x];
    Real const* AMREX_RESTRICT yoff = &voff[lo.y-vlo.y+vlen.x];

    for (int n = 0; n < ncomp; ++n) {
        for (int j = 0; j < len.y; ++j) {
            const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                fine(i,j,0,n) = crse(ic,jc,0,n)
                    + xoff[i] * slopes(ic,jc,0,n)
                    + yoff[j] * slopes(ic,jc,0,n+ncomp);
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

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                Real cmn = u(i,j,0,n);
                Real cmx = cmn;
                for     (int joff = -1; joff <= 1; ++joff) {
                    for (int ioff = -1; ioff <= 1; ++ioff) {
                        cmn = amrex::min(cmn,u(i+ioff,j+joff,0,n));
                        cmx = amrex::max(cmx,u(i+ioff,j+joff,0,n));
                    }
                }
                mm(i,j,0,n      ) = cmn - u(i,j,0,n);
                mm(i,j,0,n+ncomp) = cmx - u(i,j,0,n);
            }
        }

        compute_slopes(len, lo, hi, slopes, u, sbox, slo, shi, slen, bcr[n], n, ncomp);

        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real cen  = slopes(i,j,0,n);
                Real forw = 2.0*(u(i+1,j,0,n)-u(i  ,j,0,n));
                Real back = 2.0*(u(i  ,j,0,n)-u(i-1,j,0,n));
                Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                slopes(i,j,0,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));

                cen  = slopes(i,j,0,n+ncomp);
                forw = 2.0*(u(i,j+1,0,n)-u(i,j  ,0,n));
                back = 2.0*(u(i,j  ,0,n)-u(i,j-1,0,n));
                slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                slopes(i,j,0,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
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

    const auto clo = amrex::coarsen(lo, ratio);
    const auto slopes = ccfab.view(clo);
    const auto mm     = ccfab.view(clo,ncomp*AMREX_SPACEDIM); // min and max

    const Box& alpha_box = alphafab.box();
    const auto vlo  = amrex::lbound(alpha_box);
    const auto vlen = amrex::length(alpha_box);
    Real const* AMREX_RESTRICT xoff = &voff[lo.x-vlo.x];
    Real const* AMREX_RESTRICT yoff = &voff[lo.y-vlo.y+vlen.x];

    for (int n = 0; n < ncomp; ++n) {
        for (int j = 0; j < len.y; ++j) {
            const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                const Real dummy_fine = xoff[i]*slopes(ic,jc,0,n)
                    +                   yoff[j]*slopes(ic,jc,0,n+ncomp);

                if (dummy_fine > mm(ic,jc,0,n+ncomp) && dummy_fine != 0.0) {
                    alpha(i,j,0,n) = mm(ic,jc,0,n+ncomp) / dummy_fine;
                } else if (dummy_fine < mm(ic,jc,0,n) && dummy_fine != 0.0) {
                    alpha(i,j,0,n) = mm(ic,jc,0,n) / dummy_fine;
                } else {
                    alpha(i,j,0,n) = 1.0;
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
        for (int j = 0; j < len.y; ++j) {
            const int jj = j*ratio[1];
            for (int i = 0; i < len.x; ++i) {
                const int ii = i*ratio[0];
                Real a = 1.0;
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
        for (int j = 0; j < len.y; ++j) {
            const int jc = amrex::coarsen(j+lo.y,ratio[1]) - clo.y;
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const int ic = amrex::coarsen(i+lo.x,ratio[0]) - clo.x;
                fine(i,j,0,n) = crse(ic,jc,0,n);
            }
        }
    }
}

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int ixy  = 2;
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

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real dx0 = u(i+1,j,0,n) - u(i,j,0,n);
                Real d0x = u(i,j+1,0,n) - u(i,j,0,n);
                Real dx1 = u(i+1,j+1,0,n) - u(i,j+1,0,n);

                slope(i,j,0,n+ncomp*ix ) = rx*dx0;
                slope(i,j,0,n+ncomp*iy ) = ry*d0x;
                slope(i,j,0,n+ncomp*ixy) = rx*ry*(dx1 - dx0);
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
        for (int j = 0; j < len.y; ++j) {
            const int jc = amrex::min(amrex::coarsen(j+lo.y,ratio[1]),chi.y) - clo.y;
            const Real fy = (j+lo.y) - (jc+clo.y)*ratio[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                const int ic = amrex::min(amrex::coarsen(i+lo.x,ratio[0]),chi.x) - clo.x;
                const Real fx = (i+lo.x) - (ic+clo.x)*ratio[0];
                fine(i,j,0,n) = crse(ic,jc,0,n)
                    + fx*slope(ic,jc,0,n+ncomp*ix)
                    + fy*slope(ic,jc,0,n+ncomp*iy)
                    + fx*fy*slope(ic,jc,0,n+ncomp*ixy);
            }
        }
    }
}

}

#endif
