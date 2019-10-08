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
* \param reg
* \param rcomp
* \param flx
* \param fcomp
* \param ncomp
* \param dir
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd (Box const& bx, Array4<Real> const& reg, const int rcomp,
                 Array4<Real const> const& flx, const int fcomp, const int ncomp,
                 const int dir, Dim3 const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);

    switch (dir) {
    case 0:
    {
        const int ic = lo.x;
        const int i = ic*ratio.x;
        for (int n = 0; n < ncomp; ++n) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int joff = 0; joff < ratio.y; ++joff) {
                    const int j = ratio.y*jc + joff;
                    reg(ic,jc,0,n+rcomp) += mult*flx(i,j,0,n+fcomp);
                }
            }
        }
        break;
    }
    default:
    {
        const int jc = lo.y;
        const int j = jc*ratio.y;
        for (int n = 0; n < ncomp; ++n) {
            for (int ic = lo.x; ic <= hi.x; ++ic) {
                for (int ioff = 0; ioff < ratio.x; ++ioff) {
                    const int i = ratio.x*ic + ioff;
                    reg(ic,jc,0,n+rcomp) += mult*flx(i,j,0,n+fcomp);
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
* \param reg
* \param rcomp
* \param area
* \param flx
* \param fcomp
* \param ncomp
* \param dir
* \param ratio
* \param mult
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineareaadd (Box const& bx, Array4<Real> const& reg, const int rcomp,
                     Array4<Real const> const& area,
                     Array4<Real const> const& flx, const int fcomp, const int ncomp,
                     const int dir, Dim3 const& ratio, const Real mult) noexcept
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);

    switch (dir) {
    case 0:
    {
        const int ic = lo.x;
        const int i = ic*ratio.x;
        for (int n = 0; n < ncomp; ++n) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int joff = 0; joff < ratio.y; ++joff) {
                    const int j = ratio.y*jc + joff;
                    reg(ic,jc,0,n+rcomp) += mult*area(i,j,0)
                                                * flx(i,j,0,n+fcomp);
                }
            }
        }
        break;
    }
    default:
    {
        const int jc = lo.y;
        const int j = jc*ratio.y;
        for (int n = 0; n < ncomp; ++n) {
            for (int ic = lo.x; ic <= hi.x; ++ic) {
                for (int ioff = 0; ioff < ratio.x; ++ioff) {
                    const int i = ratio.x*ic + ioff;
                    reg(ic,jc,0,n+rcomp) += mult*area(i,j,0)
                                                * flx(i,j,0,n+fcomp);
                }
            }
        }
        break;
    }
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
        const int dir = face.coordDir();
        switch (dir) {
        case 0:
        {
            for (int n = 0; n < ncomp; ++n) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        s(i,j,0,n+scomp) += -mult*(f(i+1,j,0,n)/v(i,j,0));
                    }
                }
            }
            break;
        }
        default:
        {
            for (int n = 0; n < ncomp; ++n) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        s(i,j,0,n+scomp) += -mult*(f(i,j+1,0,n)/v(i,j,0));
                    }
                }
            }
            break;
        }
        }
    } else {
        for (int n = 0; n < ncomp; ++n) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    s(i,j,0,n+scomp) += mult*(f(i,j,0,n)/v(i,j,0));
                }
            }
        }
    }
}

}

#endif
