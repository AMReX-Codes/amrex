#ifndef AMREX_FLUXREG_2D_C_H_
#define AMREX_FLUXREG_2D_C_H_
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
* ratio is the refinement ratio. ratio must be 2.
*
* \param CrseBox
* \param reg
* \param SurfaceFluxes_Fine
* \param nFields
* \param nDOFX_X
* \param WeightsX_X
* \param LX_X
* \param ratio
*/
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd_dg( Box const&                CrseBox,
                    Array4<Real> const&       reg,
                    Array4<Real const> const& SurfaceFluxes_Fine,
                    const int                 iDimX,
                    const int                 nFields,
                    int                       nDOFX_X,
                    Real *                    WeightsX_X,
                    Real ***                  LX_X,
                    Dim3 const&               ratio ) noexcept
{
    const auto lo = amrex::lbound(CrseBox);
    const auto hi = amrex::ubound(CrseBox);

    switch( iDimX )
    {
        case 0:
        {
            const int iCrseF = lo.x;
            const int iFineF = iCrseF*ratio.x;
            int jFineV;
            int iFn;
            for( int jCrseV = lo.y; jCrseV <= hi.y; jCrseV++ ) {
                for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
                {
                    iFn = -1;
                    for( int jFn = 0; jFn < ratio.y; jFn++ ) {
                        jFineV = jCrseV*ratio.y + jFn;
                        iFn += 1;
                        for( int iNX_F = 0; iNX_F < nDOFX_X; iNX_F++ )
                        {
                            for( int iField = 0; iField < nFields; iField++ )
                            {
                                reg(iCrseF,jCrseV,0,iNX_C+iField*nDOFX_X)
                                  += WeightsX_X[iNX_F]
                                       * SurfaceFluxes_Fine
                                           (iFineF,jFineV,0,
                                            iNX_F+iField*nDOFX_X)
                                       * LX_X[iNX_C][iFn][iNX_F]
                                       * amrex::DG::FaceRatio;
                            } /* iField */
                        } /* iNX_F */
                    } /* Fine cells */
                } /* iNX_C */
            } /* Coarse cells */
            break;
        } /* case 0 */
        default:
        {
            const int jCrseF = lo.y;
            const int jFineF = jCrseF*ratio.y;
            int iFineV;
            int jFn;
            for( int iCrseV = lo.x; iCrseV <= hi.x; iCrseV++ ) {
                for( int jNX_C = 0; jNX_C < nDOFX_X; jNX_C++ )
                {
                    jFn = -1;
                    for( int iFn = 0; iFn < ratio.x; iFn++ ) {
                        iFineV = iCrseV*ratio.x + iFn;
                        jFn += 1;
                        for( int jNX_F = 0; jNX_F < nDOFX_X; jNX_F++ )
                        {
                            for( int iField = 0; iField < nFields; iField++ )
                            {
                                reg(iCrseV,jCrseF,0,jNX_C+iField*nDOFX_X)
                                  += WeightsX_X[jNX_F]
                                       * SurfaceFluxes_Fine
                                           (iFineV,jFineF,0,
                                            jNX_F+iField*nDOFX_X)
                                       * LX_X[jNX_C][jFn][jNX_F]
                                       * amrex::DG::FaceRatio;
                            } /* iField */
                        } /* jNX_F */
                    } /* Fine cells */
                } /* jNX_C */
            } /* Coarse cells */
            break;
        } /* default */
    } /* switch( iDimX ) */
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
    Array4<Real> const& G,
    Array4<Real> const& dU,
    Array4<Real const> const& dF,
    int nFields,
    Real dX1[],
    Real dX2[],
    Real dX3[],
    const Orientation face ) noexcept
{
    amrex::ignore_unused( dX3 );

    const auto CrseLo  = amrex::lbound(CrseBox);
    const auto CrseHi  = amrex::ubound(CrseBox);

    const int iDimX = face.coordDir();

    int iNX1, iNX_X1;
    int iNX2, iNX_X2;

    /* face.isLow(): face is on lower side of fine patch */

    if( face.isLow() )
    {
        switch( iDimX )
        {
            case 0:
            {
                for( int jCrse = CrseLo.y; jCrse <= CrseHi.y; ++jCrse ) {
                for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
                    for( int iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                    {
                        iNX1   = amrex::DG::NodeNumberTableX[0][iNX];
                        iNX_X1 = amrex::DG::NodeNumberTableX_X1[iNX];
                        for( int iField = 0; iField < nFields; iField++ )
                        {
                            dU(iCrse,jCrse,0,iNX+iField*amrex::DG::nDOFX)
                              += -dF(iCrse+1,jCrse,0,
                                     iNX_X1+iField*amrex::DG::nDOFX_X1)
                                   * amrex::DG::LX_X1_Up[iNX1]
                                   / ( amrex::DG::WeightsX_q[iNX]
                                         * G(iCrse,jCrse,0,
                                             (amrex::DG::iGF_SqrtGm-1)
                                                *amrex::DG::nDOFX+iNX)
                                         * dX1[iCrse] );
                        } /* iField */
                    } /* iNX */
                }} /* Coarse cells */
                break;
            } /* case 0 */
            default:
            {
                for( int jCrse = CrseLo.y; jCrse <= CrseHi.y; ++jCrse ) {
                for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
                    for( int iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                    {
                        iNX2   = amrex::DG::NodeNumberTableX[1][iNX];
                        iNX_X2 = amrex::DG::NodeNumberTableX_X2[iNX];
                        for( int iField = 0; iField < nFields; iField++ )
                        {
                            dU(iCrse,jCrse,0,iNX+iField*amrex::DG::nDOFX)
                              += -dF(iCrse,jCrse+1,0,
                                     iNX_X2+iField*amrex::DG::nDOFX_X2)
                                   * amrex::DG::LX_X2_Up[iNX2]
                                   / ( amrex::DG::WeightsX_q[iNX]
                                         * G(iCrse,jCrse,0,
                                             (amrex::DG::iGF_SqrtGm-1)
                                                *amrex::DG::nDOFX+iNX)
                                         * dX2[jCrse] );
                        } /* iField */
                    } /* iNX */
                }} /* Coarse cells */
                break;
            } /* default */
        } /* switch( iDimX ) */
    }
    else /* face is on upper side of fine patch */
    {
        switch( iDimX )
        {
            case 0:
            {
                for( int jCrse = CrseLo.y; jCrse <= CrseHi.y; ++jCrse ) {
                for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
                    for( int iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                    {
                        iNX1   = amrex::DG::NodeNumberTableX[0][iNX];
                        iNX_X1 = amrex::DG::NodeNumberTableX_X1[iNX];
                        for( int iField = 0; iField < nFields; iField++ )
                        {
                            dU(iCrse,jCrse,0,iNX+iField*amrex::DG::nDOFX)
                              += dF(iCrse,jCrse,0,
                                    iNX_X1+iField*amrex::DG::nDOFX_X1)
                                   * amrex::DG::LX_X1_Dn[iNX1]
                                   / ( amrex::DG::WeightsX_q[iNX]
                                         * G(iCrse,jCrse,0,
                                             (amrex::DG::iGF_SqrtGm-1)
                                                *amrex::DG::nDOFX+iNX)
                                         * dX1[iCrse] );
                        } /* iField */
                    } /* iNX */
                }} /* Coarse cells */
                break;
            } /* case 0 */
            default:
            {
                for( int jCrse = CrseLo.y; jCrse <= CrseHi.y; ++jCrse ) {
                for( int iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse ) {
                    for( int iNX = 0; iNX < amrex::DG::nDOFX; iNX++ )
                    {
                        iNX2   = amrex::DG::NodeNumberTableX[1][iNX];
                        iNX_X2 = amrex::DG::NodeNumberTableX_X2[iNX];
                        for( int iField = 0; iField < nFields; iField++ )
                        {
                            dU(iCrse,jCrse,0,iNX+iField*amrex::DG::nDOFX)
                              += dF(iCrse,jCrse,0,
                                     iNX_X2+iField*amrex::DG::nDOFX_X2)
                                   * amrex::DG::LX_X2_Dn[iNX2]
                                   / ( amrex::DG::WeightsX_q[iNX]
                                         * G(iCrse,jCrse,0,
                                             (amrex::DG::iGF_SqrtGm-1)
                                                *amrex::DG::nDOFX+iNX)
                                         * dX2[jCrse] );
                        } /* iField */
                    } /* iNX */
                }} /* Coarse elements */
                break;
            } /* case default */
        } /* switch( iDimX ) */
    }
} /* END fluxreg_reflux_dg */


}

#endif
