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
* \brief Add fine grid flux to flux register (DG).  Flux array is a fine grid
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
fluxreg_fineadd_dg( Box const& bx,
                    Array4<Real> const& reg,
                    Array4<Real const> const& flx,
                    const int /*dir*/,
                    const int nFields,
                    const int nDOFX_X,
                    const Real WeightsX_X[],
                    const Real LX_X[],
                    Dim3 const& ratio,
                    const Real mult ) noexcept
{
    const Real Two = 2.0;

    const int nFine = pow(2,AMREX_SPACEDIM-1);

    /* Convert LX_X into 3D array */

    Real *** LX_XX = new Real ** [nDOFX_X];
    for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
    {
        LX_XX[iNX_C] = new Real * [nFine];

        for( int iFine = 0; iFine < nFine; iFine++ )
        {
          LX_XX[iNX_C][iFine] = new Real [nDOFX_X];
        }
    }

    int k = -1;
    for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ ) {
    for( int iFine = 0; iFine < nFine  ; iFine++ ) {
    for( int iNX_F = 0; iNX_F < nDOFX_X; iNX_F++ ) {
        k += 1;
        LX_XX[iNX_C][iFine][iNX_F] = LX_X[k];
    }}}

    const auto lo   = amrex::lbound(bx);
    const int iCrse = lo.x;
    const int iFine = iCrse*ratio.x;

    /* Need to differentiate iFine being on lower or higher side of iCrse */

    for( int iField = 0; iField < nFields; iField++ )
    {
        for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
        {
            for( int iNX_F = 0; iNX_F < nDOFX_X; iNX_F++ )
            {
                reg(iCrse,0,0,iNX_C+iField*nDOFX_X)
                  += mult * Two * WeightsX_X[iNX_F]
                       * flx(iFine,0,0,iNX_F+iField*nDOFX_X)
                       * LX_XX[iNX_C][0][iNX_F];
            } /* iNX_F */
        } /* iNX_C */
    } /* iField */

    for( int iNX_C = 0; iNX_C < nDOFX_X; iNX_C++ )
    {
        for( int iFn = 0; iFn < nFine; iFn++ )
        {
            delete LX_XX[iNX_C][iFn];
        }
        delete LX_XX[iNX_C];
    }
    delete LX_XX;

//    for( int n = 0; n < ncomp; ++n )
//    {
//        reg(iCrse,0,0,n+rcomp) += mult * flx(iFine,0,0,n+fcomp);
//    }
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


AMREX_GPU_HOST_DEVICE inline void
fluxreg_reflux_dg
  ( Box const& CrseBox,
    Array4<Real> const& u,
    Array4<Real const> const& f,
    const Orientation face,
    int nFields,
    int nDOFX,
    int nNodesX[],
    Real WeightsX_q[],
    Real dX1[], Real dX2[], Real dX3[],
    Real LX_X1_Dn[], Real LX_X1_Up[],
    Real LX_X2_Dn[], Real LX_X2_Up[],
    Real LX_X3_Dn[], Real LX_X3_Up[] ) noexcept
{

    const auto CrseLo  = amrex::lbound(CrseBox);
    const auto CrseHi  = amrex::ubound(CrseBox);

    const int nDOFX_X = nNodesX[1] * nNodesX[2];

    int iNX, iNX1, iNX2, iNX3, iNX_X;
    int iCrse, iField;

    int ** NodeNumberTableX = new int * [3];
    for( int iDimX = 0; iDimX < 3; iDimX++ )
    {
        NodeNumberTableX[iDimX] = new int [nDOFX];
    }

    int * NodeNumberTableX_X = new int [nDOFX];

    int *** NodeNumberTableX3D = new int ** [nNodesX[0]];
    for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ )
    {
        NodeNumberTableX3D[iNX1] = new int * [nNodesX[1]];
        for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ )
        {
            NodeNumberTableX3D[iNX1][iNX2] = new int [nNodesX[2]];
        }
    }

    iNX = -1;
    for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      iNX += 1;
      NodeNumberTableX[0][iNX] = iNX1;
      NodeNumberTableX[1][iNX] = iNX2;
      NodeNumberTableX[2][iNX] = iNX3;
      NodeNumberTableX3D[iNX1][iNX2][iNX3] = iNX;
    }}}

    for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      iNX = NodeNumberTableX3D[iNX1][iNX2][iNX3];
      NodeNumberTableX_X[iNX] = iNX / nNodesX[0];
    }}}

    /* Populate destination MultiFab */
    if( face.isLow() )
    {
        for( iCrse = CrseLo.x; iCrse <= CrseHi.x; ++iCrse )
        {
            for( iField = 0; iField < nFields; iField++ )
            {
                for( iNX = 0; iNX < nDOFX; iNX++ )
                {
                    iNX1  = NodeNumberTableX[0][iNX];
                    iNX_X = NodeNumberTableX_X[iNX];

                    u(iCrse,0,0,iNX+iField*nDOFX)
                      += -f(iCrse+1,0,0,iNX_X+iField*nDOFX_X)
                           * LX_X1_Up[iNX1] / ( WeightsX_q[iNX] * dX1[iCrse] );
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
                for( iNX = 0; iNX < nDOFX; iNX++ )
                {
                    iNX1  = NodeNumberTableX[0][iNX];
                    iNX_X = NodeNumberTableX_X[iNX];

                    u(iCrse,0,0,iNX+iField*nDOFX)
                      += f(iCrse,0,0,iNX_X+iField*nDOFX_X)
                           * LX_X1_Dn[iNX1] / ( WeightsX_q[iNX] * dX1[iCrse] );
                }
            }
        }
    }

    for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ )
    {
        for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ )
        {
            delete NodeNumberTableX3D[iNX1][iNX2];
        }
        delete NodeNumberTableX3D[iNX1];
    }
    delete NodeNumberTableX3D;

    delete NodeNumberTableX_X;
    for( int iDimX = 0; iDimX < 3; iDimX++ )
    {
        delete NodeNumberTableX[iDimX];
    }
    delete NodeNumberTableX;
} /* END fluxreg_reflux_dg */

}

#endif
