#ifndef AMREX_FLUXREG_2D_C_H_
#define AMREX_FLUXREG_2D_C_H_
#include <AMReX_Config.H>

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
* \brief Add fine grid flux to flux register (DG).
* bx is a box defined on the coarse level.
* reg is a coarse grid edge based object.
* SurfaceFluxes_Fine is a fine grid edge based object. It is assumed that the
* coarsened flux region contains the register region.
* iDimX is the spatial dimension in which the fluxes are defined.
* nFields is the number of fields in the grid based objects.
* nDOFX_X is the number of degrees of freedom per field on the interface in the
* iDimX dimension.
* WeightsX_X are the Gaussian quadrature weights on the interface in the iDimX
* dimension. These Gaussian points are defined on the interval [-0.5,+0.5].
* LX_X are the Lagrange polynomials evaluated at the quadrature points in the
* AMREX_SPACDIM-1 dimensions on the interface in the iDimX dimension.
* ratio is the refinement ratio. ratio must be 2.
*
* \param bx
* \param reg
* \param SurfaceFluxes_Fine
* \param iDimX
* \param nFields
* \param nDOFX_X
* \param WeightsX_X
* \param LX_X
* \param ratio
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd_dg( Box const&                bx,
                    Array4<Real> const&       reg,
                    Array4<Real const> const& SurfaceFluxes_Fine,
                    const int                 iDimX,
                    const int                 nFields,
                    int                       nDOFX_X,
                    Real *                    WeightsX_X,
                    Real ***                  LX_X,
                    Dim3 const&               ratio ) noexcept
{


#ifndef AMREX_USE_GPU
    std::cout<<"\n\nfluxreg_fineadd_dg not implemented for AMREX_SPACEDIM=2\n\n";
#endif

    amrex::ignore_unused
      ( bx, reg, SurfaceFluxes_Fine, iDimX, nFields, nDOFX_X, WeightsX_X, LX_X, ratio );

} /* END fluxreg_fineadd_dg */


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

/**
* \brief Replace coarse flux with fine flux (DG).
* CrseBox is a box defined on the coarse level.
* dU is the increment for each DOF within an element.
* dF is a coarse grid edge based object.
* nFields is the number of fields in the grid based objects.
* dX1 is the width of the coarse element in the X1-dimension.
* dX2 is the width of the coarse element in the X2-dimension.
* dX3 is the width of the coarse element in the X3-dimension.
* face is ...
*
* \param CrseBox
* \param dU
* \param dF
* \param nFields
* \param dX1
* \param dX2
* \param dX3
* \param face
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_reflux_dg
  ( Box const& CrseBox,
    Array4<Real> const& dU,
    Array4<Real const> const& dF,
    int nFields,
    Real dX1[],
    Real dX2[],
    Real dX3[],
    const Orientation face ) noexcept
{

#ifndef AMREX_USE_GPU
    std::cout<<"\n\nfluxreg_reflux_dg not implemented for AMREX_SPACEDIM=2\n\n";
#endif

    amrex::ignore_unused
      ( CrseBox, dU, dF, nFields, dX1, dX2, dX3, face );

} /* END fluxreg_reflux_dg */


}

#endif
