#ifndef AMREX_FLUXREG_2D_C_H_
#define AMREX_FLUXREG_2D_C_H_

#include <AMReX_FArrayBox.H>

namespace amrex {


/**
* \brief Add fine grid flux to flux register.  Flux array is a fine grid
* edge based object, Register is a coarse grid edge based object.
* It is assumed that the coarsened flux region contains the register
* region.
*
* \param bx
* \param regfab
* \param rcomp
* \param flxfab
* \param fcomp
* \param ncomp
* \param dir
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd (Box const& bx, FArrayBox& regfab, const int rcomp,
                 FArrayBox const& flxfab, const int fcomp, const int ncomp,
                 const int dir, IntVect const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto reg = regfab.view(lo,rcomp);

    const auto flo = amrex::lbound(flxfab.box());
    const auto flx = flxfab.view(flo,fcomp);

    switch (dir) {
    case 0:
    {
        const int ic = lo.x;
        const int i = ic*ratio[0];
        for (int n = 0; n < ncomp; ++n) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int joff = 0; joff < ratio[1]; ++joff) {
                    const int j = ratio[1]*jc + joff;
                    reg(ic-lo.x,jc-lo.y,0,n) += mult
                        *flx(i-flo.x,j-flo.y,0,n);
                }
            }
        }
        break;
    }
    default:
    {
        const int jc = lo.y;
        const int j = jc*ratio[1];
        for (int n = 0; n < ncomp; ++n) {
            for (int ic = lo.x; ic <= hi.x; ++ic) {
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    const int i = ratio[0]*ic + ioff;
                    reg(ic-lo.x,jc-lo.y,0,n) += mult
                        *flx(i-flo.x,j-flo.y,0,n);
                }
            }
        }
    }
    }
}


/**
* \brief Add fine grid flux times area to flux register.  Flux array is a fine grid
* edge based object, Register is a coarse grid edge based object.
* It is assumed that the coarsened flux region contains the register
* region.
*
* \param bx
* \param regfab
* \param rcomp
* \param areafab
* \param flxfab
* \param fcomp
* \param ncomp
* \param dir
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineareaadd (Box const& bx, FArrayBox& regfab, const int rcomp,
                     FArrayBox const& areafab,
                     FArrayBox const& flxfab, const int fcomp, const int ncomp,
                     const int dir, IntVect const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto reg = regfab.view(lo,rcomp);

    const auto flo = amrex::lbound(flxfab.box());
    const auto flx = flxfab.view(flo,fcomp);
    const auto area = areafab.view(flo);

    switch (dir) {
    case 0:
    {
        const int ic = lo.x;
        const int i = ic*ratio[0];
        for (int n = 0; n < ncomp; ++n) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int joff = 0; joff < ratio[1]; ++joff) {
                    const int j = ratio[1]*jc + joff;
                    reg(ic-lo.x,jc-lo.y,0,n) += mult
                        *area(i-flo.x,j-flo.y,0)
                        * flx(i-flo.x,j-flo.y,0,n);
                }
            }
        }
        break;
    }
    default:
    {
        const int jc = lo.y;
        const int j = jc*ratio[1];
        for (int n = 0; n < ncomp; ++n) {
            for (int ic = lo.x; ic <= hi.x; ++ic) {
                for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                    const int i = ratio[0]*ic + ioff;
                    reg(ic-lo.x,jc-lo.y,0,n) += mult
                        *area(i-flo.x,j-flo.y,0)
                        * flx(i-flo.x,j-flo.y,0,n);
                }
            }
        }
        break;
    }
    }
}

AMREX_GPU_HOST_DEVICE inline void
fluxreg_reflux (Box const& bx, FArrayBox& sfab, const int scomp,
                FArrayBox const& ffab, FArrayBox const& vfab,
                const int ncomp, const Real mult, const Orientation face) noexcept
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto s = sfab.view(lo,scomp);
    const auto f = ffab.view(lo);
    const auto v = vfab.view(lo);

    if (face.isLow()) {
        const int dir = face.coordDir();
        switch (dir) {
        case 0:
        {
            for (int n = 0; n < ncomp; ++n) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        s(i,j,0,n) += -mult*(f(i+1,j,0,n)/v(i,j,0));
                    }
                }
            }
            break;
        }
        default:
        {
            for (int n = 0; n < ncomp; ++n) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        s(i,j,0,n) += -mult*(f(i,j+1,0,n)/v(i,j,0));
                    }
                }
            }
            break;
        }
        }
    } else {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    s(i,j,0,n) += mult*(f(i,j,0,n)/v(i,j,0));
                }
            }
        }
    }
}

}

#endif
