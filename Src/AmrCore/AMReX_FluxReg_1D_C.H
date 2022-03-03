#ifndef AMREX_FLUXREG_1D_C_H_
#define AMREX_FLUXREG_1D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_DG.H>

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
* \brief Add fine grid flux to flux register (DG).
* CrseBx is a box defined on the coarse level.
* reg is a coarse grid edge based object.
* SurfaceFluxes_Fine is a fine grid edge based object. It is assumed that the
* coarsened flux region contains the register region.
* nFields is the number of fields in the grid based objects.
* nDOFX_X is the number of degrees of freedom per field on the interface in the
* iDimX dimension.
* WeightsX_X are the Gaussian quadrature weights on the interface in the iDimX
* dimension. These Gaussian points are defined on the interval [-0.5,+0.5].
* LX_X are the Lagrange polynomials evaluated at the quadrature points in the
* AMREX_SPACDIM-1 dimensions on the interface in the iDimX dimension.
* ratio is the refinement ratio. ratio must be 2.
*
* \param CrseBx
* \param reg
* \param SurfaceFluxes_Fine
* \param nFields
* \param nDOFX_X
* \param WeightsX_X
* \param LX_X
* \param ratio
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd_dg( Box const&                CrseBx,
                    Array4<Real> const&       reg,
                    Array4<Real const> const& SurfaceFluxes_Fine,
                    const int                 /*iDimX*/,
                    const int                 nFields,
                    int                       nDOFX_X,
                    Real *                    WeightsX_X,
                    Real ***                  LX_X,
                    Dim3 const&               ratio ) noexcept
{

    const auto lo   = amrex::lbound(CrseBx);
    const int iCrse = lo.x;
    const int iFine = iCrse*ratio.x;

    for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
    {
        for( int iFn = 0; iFn < amrex::DG::nFineF; iFn++ )
        {
            for( int iField = 0; iField < nFields; iField++ )
            {
                for( int iNX_F = 0; iNX_F < nDOFX_X; iNX_F++ )
                {
                    reg(iCrse,0,0,iNX_C+iField*nDOFX_X)
                      += WeightsX_X[iNX_F]
                           * SurfaceFluxes_Fine(iFine,0,0,
                                  iNX_F+iField*nDOFX_X+nFields*nDOFX_X*iFn)
                           * LX_X[iNX_C][iFn][iNX_F]
                           * amrex::DG::FaceRatio;
                } /* iNX_F */
            } /* iField */
        } /* iFn */
    } /* iNX_C */

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

    amrex::ignore_unused( dX2, dX3 );

    const auto CrseLo  = amrex::lbound(CrseBox);
    const auto CrseHi  = amrex::ubound(CrseBox);

    int iNX, iNX1, iNX_X1;
    int iCrse, iField;

    /* face.isLow(): face is on lower side of fine patch */

    if( face.isLow() )
    {
        for( iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse )
        {
            for( iField = 0; iField < nFields; iField++ )
            {
                for( iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                {
                    iNX1   = amrex::DG::NodeNumberTableX[0][iNX];
                    iNX_X1 = amrex::DG::NodeNumberTableX_X1[iNX];

                    dU(iCrse,0,0,iNX+iField*amrex::DG::nDOFX)
                      += -dF(iCrse+1,0,0,iNX_X1+iField*amrex::DG::nDOFX_X1)
                           * amrex::DG::LX_X1_Up[iNX1]
                           / ( amrex::DG::WeightsX_q[iNX] * dX1[iCrse] );
                }
            }
        }
    }
    else
    {
        for( iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse )
        {
            for( iField = 0; iField < nFields; iField++ )
            {
                for( iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                {
                    iNX1   = amrex::DG::NodeNumberTableX[0][iNX];
                    iNX_X1 = amrex::DG::NodeNumberTableX_X1[iNX];

                    dU(iCrse,0,0,iNX+iField*amrex::DG::nDOFX)
                      += dF(iCrse,0,0,iNX_X1+iField*amrex::DG::nDOFX_X1)
                           * amrex::DG::LX_X1_Dn[iNX1]
                           / ( amrex::DG::WeightsX_q[iNX] * dX1[iCrse] );
                }
            }
        }
    }

} /* END fluxreg_reflux_dg */

}

#endif
