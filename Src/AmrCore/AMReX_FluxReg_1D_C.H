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
* \param reg
* \param rcomp
* \param flx
* \param fcomp
* \param ncomp
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd (Box const& bx, Array4<Real> const& reg, const int rcomp,
                 Array4<Real const> const& flx, const int fcomp, const int ncomp,
                 const int /*dir*/, Dim3 const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const int ic = lo.x;
    const int i = ic*ratio.x;
    for (int n = 0; n < ncomp; ++n) {
        reg(ic,0,0,n+rcomp) += mult * flx(i,0,0,n+fcomp);
    }
}


/**
* \brief Add fine grid flux times area to flux register.  Flux array is a fine grid
* edge based object, Register is a coarse grid edge based object.
* It is assumed that the coarsened flux region contains the register
* region.
*
* \param bx
* \param reg
* \param rcomp
* \param area
* \param flx
* \param fcomp
* \param ncomp
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineareaadd (Box const& bx, Array4<Real> const& reg, const int rcomp,
                     Array4<Real const> const& area,
                     Array4<Real const> const& flx, const int fcomp, const int ncomp,
                     const int /*dir*/, Dim3 const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const int ic = lo.x;
    const int i = ic*ratio.x;
    for (int n = 0; n < ncomp; ++n) {
        reg(ic,0,0,n+rcomp) += mult * area(i,0,0) * flx(i,0,0,n+fcomp);
    }
}

AMREX_GPU_HOST_DEVICE inline void
fluxreg_reflux (Box const& bx, Array4<Real> const& s, const int scomp,
                Array4<Real const> const& f, Array4<Real const> const& v,
                const int ncomp, const Real mult, const Orientation face) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    if (face.isLow()) {
        for (int n = 0; n < ncomp; ++n) {
            for (int i = lo.x; i <= hi.x; ++i) {
                s(i,0,0,n+scomp) += -mult*f(i+1,0,0,n)/v(i,0,0);
            }
        }
    } else {
        for (int n = 0; n < ncomp; ++n) {
            for (int i = lo.x; i <= hi.x; ++i) {
                s(i,0,0,n+scomp) += mult*f(i,0,0,n)/v(i,0,0);
            }
        }
    }
}

}

#endif
