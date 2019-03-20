#ifndef AMREX_FLUXREG_1D_C_H_
#define AMREX_FLUXREG_1D_C_H_

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
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd (Box const& bx, FArrayBox& regfab, const int rcomp,
                 FArrayBox const& flxfab, const int fcomp, const int ncomp,
                 const int /*dir*/, IntVect const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto reg = regfab.view(lo,rcomp);

    const auto flo = amrex::lbound(flxfab.box());
    const auto flx = flxfab.view(flo,fcomp);

    const int ic = lo.x;
    const int i = ic*ratio[0];
    for (int n = 0; n < ncomp; ++n) {
        reg(ic-lo.x,0,0,n) += mult * flx(i-flo.x,0,0,n);
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
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineareaadd (Box const& bx, FArrayBox& regfab, const int rcomp,
                     FArrayBox const& areafab,
                     FArrayBox const& flxfab, const int fcomp, const int ncomp,
                     const int /*dir*/, IntVect const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto reg = regfab.view(lo,rcomp);

    const auto flo = amrex::lbound(flxfab.box());
    const auto flx = flxfab.view(flo,fcomp);
    const auto area = areafab.view(flo);

    const int ic = lo.x;
    const int i = ic*ratio[0];
    for (int n = 0; n < ncomp; ++n) {
        reg(ic-lo.x,0,0,n) += mult * area(i-flo.x,0,0) * flx(i-flo.x,0,0,n);
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
        for (int n = 0; n < ncomp; ++n) {
            for (int i = 0; i < len.x; ++i) {
                s(i,0,0,n) += -mult*f(i+1,0,0,n)/v(i,0,0);
            }
        }
    } else {
        for (int n = 0; n < ncomp; ++n) {
            for (int i = 0; i < len.x; ++i) {
                s(i,0,0,n) += mult*f(i,0,0,n)/v(i,0,0);
            }
        }
    }
}

}

#endif
