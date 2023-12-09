#ifndef AMREX_FLUXREG_1D_C_H_
#define AMREX_FLUXREG_1D_C_H_
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
* CrseBox is a box defined on the coarse level.
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
* FaceRatio is the ratio of a fine face to a coarse face
* ratio is the refinement ratio. ratio must be 2.
*
* \param CrseBox
* \param reg
* \param SurfaceFluxes_Fine
* \param nFields
* \param nDOFX_X
* \param WeightsX_X
* \param LX_X
* \param FaceRatio
* \param ratio
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd_dg( Box const&                CrseBox,
                    Array4<Real> const&       reg,
                    Array4<Real const> const& SurfaceFluxes_Fine,
                    const int                 /*iDimX*/,
                    const int                 nFields,
                    int                       nDOFX_X,
                    const Real *              WeightsX_X,
                    Array4<Real const>        LX_X,
                    Real                      FaceRatio,
                    Dim3 const&               ratio ) noexcept
{
    const auto lo    = amrex::lbound(CrseBox);
    const int iCrseF = lo.x;
    const int iFineF = iCrseF*ratio.x;

    for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
    {
        for( int iNX_F = 0; iNX_F < nDOFX_X; iNX_F++ )
        {
            for( int iField = 0; iField < nFields; iField++ )
            {
                reg(iCrseF,0,0,iNX_C+iField*nDOFX_X)
                  += WeightsX_X[iNX_F]
                       * SurfaceFluxes_Fine(iFineF,0,0,iNX_F+iField*nDOFX_X)
                       * LX_X(0,iNX_C,0,iNX_F)
                       * FaceRatio;
            } /* iField */
        } /* iNX_F */
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
* G is the array of geometry fields, in particular the square root of the
*   spatial three-metric determinant.
* dU is the increment for each degree of freedom of each fluid field
*   within an element.
* dF is a coarse grid edge based object that holds the flux correction.
* nDOFX is the number of degrees of freedom per element, per field
* nDOFX_X1 is the number of degrees of freedom per element X1-face, per field.
* nDOFX_X2 is the number of degrees of freedom per element X2-face, per field.
* nDOFX_X3 is the number of degrees of freedom per element X3-face, per field.
* nFields is the number of fields in the grid based objects.
* iGF_SqrtGm is index of G corresponding the the square root of the spatial
*   three-metric determinant.
* NodeNumberTableX_X1 is the table that maps degrees of freedom within an
*   element to degrees of freedom on the X1 interface.
* NodeNumberTableX_X2 is the table that maps degrees of freedom within an
*   element to degrees of freedom on the X2 interface.
* NodeNumberTableX_X3 is the table that maps degrees of freedom within an
*   element to degrees of freedom on the X3 interface.
* WeightsX_q are the Gauss--Legendre quadrature weights within an element.
* LX_X1_Up are the Lagrange polynomials interpolated to the upper side of
*   an element in the X1 dimension.
* LX_X1_Dn are the Lagrange polynomials interpolated to the lower side of
*   an element in the X1 dimension.
* LX_X2_Up are the Lagrange polynomials interpolated to the upper side of
*   an element in the X2 dimension.
* LX_X2_Dn are the Lagrange polynomials interpolated to the lower side of
*   an element in the X2 dimension.
* LX_X3_Up are the Lagrange polynomials interpolated to the upper side of
*   an element in the X3 dimension.
* LX_X3_Dn are the Lagrange polynomials interpolated to the lower side of
*   an element in the X3 dimension.
* dX1 is the mesh width in the X1 dimension.
* dX2 is the mesh width in the X2 dimension.
* dX3 is the mesh width in the X3 dimension.
* face is ...
*
* \param CrseBox
* \param G
* \param dU
* \param dF
* \param nDOFX
* \param nDOFX_X1
* \param nDOFX_X2
* \param nDOFX_X3
* \param nFields
* \param iGF_SqrtGm
* \param NodeNumberTableX_X1
* \param NodeNumberTableX_X2
* \param NodeNumberTableX_X3
* \param WeightsX_q
* \param LX_X1_Up
* \param LX_X1_Dn
* \param LX_X2_Up
* \param LX_X2_Dn
* \param LX_X3_Up
* \param LX_X3_Dn
* \param dX1
* \param dX2
* \param dX3
* \param face
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_reflux_dg
  ( Box                const& CrseBox,
    Array4<Real>       const& G,
    Array4<Real>       const& dU,
    Array4<Real const> const& dF,
    int                       nDOFX,
    int                       nDOFX_X1,
    int                       /*nDOFX_X2*/,
    int                       /*nDOFX_X3*/,
    int                       nFields,
    int                       iGF_SqrtGm,
    Array4<int>               NodeNumberTableX_X1,
    Array4<int>               /*NodeNumberTableX_X2*/,
    Array4<int>               /*NodeNumberTableX_X3*/,
    Array4<Real>              WeightsX_q,
    Array4<Real>              LX_X1_Up,
    Array4<Real>              LX_X1_Dn,
    Array4<Real>              /*LX_X2_Up*/,
    Array4<Real>              /*LX_X2_Dn*/,
    Array4<Real>              /*LX_X3_Up*/,
    Array4<Real>              /*LX_X3_Dn*/,
    Real                      dX1,
    Real                      /*dX2*/,
    Real                      /*dX3*/,
    const Orientation         face ) noexcept
{

    const auto CrseLo  = amrex::lbound(CrseBox);
    const auto CrseHi  = amrex::ubound(CrseBox);

    int iNX_X1;

    /* face.isLow(): face is on lower side of fine patch */

    if( face.isLow() )
    {
        for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
            for( int iNX = 0; iNX < nDOFX; iNX++ ) {
                iNX_X1 = NodeNumberTableX_X1(iNX,0,0,0);
                for( int iField = 0; iField < nFields; iField++ )
                {
                    dU(iCrse,0,0,iNX+iField*nDOFX)
                      += -dF(iCrse+1,0,0,
                             iNX_X1+iField*nDOFX_X1)
                           * LX_X1_Up(iNX_X1,iNX,0,0)
                           / ( WeightsX_q(iNX,0,0,0)
                                 * G(iCrse,0,0,
                                     (iGF_SqrtGm-1)*nDOFX+iNX)
                                 * dX1 );
                } /* iField */
            } /* iNX */
        } /* Coarse cells */
    }
    else
    {
        for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
            for( int iNX = 0; iNX < nDOFX; iNX++ )
            {
                iNX_X1 = NodeNumberTableX_X1(iNX,0,0,0);
                for( int iField = 0; iField < nFields; iField++ )
                {
                    dU(iCrse,0,0,iNX+iField*nDOFX)
                      += dF(iCrse,0,0,
                            iNX_X1+iField*nDOFX_X1)
                          * LX_X1_Dn(iNX_X1,iNX,0,0)
                          / ( WeightsX_q(iNX,0,0,0)
                                * G(iCrse,0,0,
                                    (iGF_SqrtGm-1)*nDOFX+iNX)
                                * dX1 );
                } /* iField */
            } /* iNX */
        } /* Coarse cells */
    }
} /* END fluxreg_reflux_dg */

}

#endif
