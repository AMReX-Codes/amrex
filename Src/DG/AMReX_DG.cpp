#include <AMReX_DG.H>

#include <cstddef> /* For NULL */
#include <iostream> /* for std::cout/std::endl */
#include <cmath> /* for pow */

#include <AMReX_REAL.H>

namespace amrex
{
namespace DG
{

int nFineV;
int nFineF;

Real VolumeRatio;
Real FaceRatio;

int nNodesX[3];
int nDOFX;
int nDOFX_X1;
int nDOFX_X2;
int nDOFX_X3;

Real * WeightsX_X1 = NULL;
Real * WeightsX_X2 = NULL;
Real * WeightsX_X3 = NULL;
Real * WeightsX_q  = NULL;

Real *** ProjectionMatrix   = NULL;
Real *** ProjectionMatrix_T = NULL;

Real *** LX_X1 = NULL;
Real *** LX_X2 = NULL;
Real *** LX_X3 = NULL;

Real * LX_X1_Up = NULL;
Real * LX_X1_Dn = NULL;
Real * LX_X2_Up = NULL;
Real * LX_X2_Dn = NULL;
Real * LX_X3_Up = NULL;
Real * LX_X3_Dn = NULL;

int **  NodeNumberTableX    = NULL;
int *** NodeNumberTableX3D  = NULL;
int *   NodeNumberTableX_X1 = NULL;
int *   NodeNumberTableX_X2 = NULL;
int *   NodeNumberTableX_X3 = NULL;

int iGF_SqrtGm;

void InitializeMeshRefinement_DG
       ( int N[], Real ProjMatrix[],
         Real WeightsX1[], Real WeightsX2[], Real WeightsX3[],
         Real LX_X1_Refined_Packed[],
         Real LX_X2_Refined_Packed[],
         Real LX_X3_Refined_Packed[],
         Real LX_X1_Up_1D[], Real LX_X1_Dn_1D[],
         Real LX_X2_Up_1D[], Real LX_X2_Dn_1D[],
         Real LX_X3_Up_1D[], Real LX_X3_Dn_1D[], int iGF_SqtGm )
{
    int k;

    nFineV = (int)pow( 2, AMREX_SPACEDIM );
    nFineF = (int)pow( 2, AMREX_SPACEDIM-1 );

    VolumeRatio = One / (Real)nFineV;
    FaceRatio   = One / (Real)nFineF;

    nNodesX[0] = N[0];
    nNodesX[1] = N[1];
    nNodesX[2] = N[2];

    nDOFX = nNodesX[0] * nNodesX[1] * nNodesX[2];

    nDOFX_X1 = nNodesX[1] * nNodesX[2];
    nDOFX_X2 = nNodesX[0] * nNodesX[2];
    nDOFX_X3 = nNodesX[0] * nNodesX[1];

    AllocateArray( nDOFX_X1, WeightsX_X1 );
    AllocateArray( nDOFX_X2, WeightsX_X2 );
    AllocateArray( nDOFX_X3, WeightsX_X3 );
    AllocateArray( nDOFX, WeightsX_q );

    AllocateArray( nFineV, nDOFX, nDOFX, ProjectionMatrix   );
    AllocateArray( nFineV, nDOFX, nDOFX, ProjectionMatrix_T );

    AllocateArray( nDOFX_X1, nFineF, nDOFX_X1, LX_X1 );
    AllocateArray( nDOFX_X2, nFineF, nDOFX_X2, LX_X2 );
    AllocateArray( nDOFX_X3, nFineF, nDOFX_X3, LX_X3 );

    AllocateArray( nDOFX_X1, LX_X1_Up );
    AllocateArray( nDOFX_X1, LX_X1_Dn );
    AllocateArray( nDOFX_X2, LX_X2_Up );
    AllocateArray( nDOFX_X2, LX_X2_Dn );
    AllocateArray( nDOFX_X3, LX_X3_Up );
    AllocateArray( nDOFX_X3, LX_X3_Dn );

    AllocateArray( 3, nDOFX, NodeNumberTableX );
    AllocateArray( nNodesX[0], nNodesX[1], nNodesX[2], NodeNumberTableX3D );
    AllocateArray( nDOFX, NodeNumberTableX_X1 );
    AllocateArray( nDOFX, NodeNumberTableX_X2 );
    AllocateArray( nDOFX, NodeNumberTableX_X3 );

    k = -1;
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
      k += 1;
      WeightsX_X1[k] = WeightsX2[iNX2] * WeightsX3[iNX3];
    }}

    k = -1;
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      k += 1;
      WeightsX_X2[k] = WeightsX1[iNX1] * WeightsX3[iNX3];
    }}

    k = -1;
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      k += 1;
      WeightsX_X3[k] = WeightsX1[iNX1] * WeightsX2[iNX2];
    }}

    k = -1;
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      k += 1;
      WeightsX_q[k] = WeightsX1[iNX1] * WeightsX2[iNX2] * WeightsX3[iNX3];
    }}}

    k = -1;
    for( int iFine = 0; iFine < nFineV; iFine++ ) {
    for( int jNX   = 0; jNX   < nDOFX; jNX++   ) {
    for( int iNX   = 0; iNX   < nDOFX; iNX++   ) {
        k += 1;
        ProjectionMatrix  [iFine][iNX][jNX] = ProjMatrix[k];
        ProjectionMatrix_T[iFine][jNX][iNX] = ProjMatrix[k];
    }}}

    k = -1;
    for( int iNX_C = 0; iNX_C < nDOFX_X1; iNX_C++ ) {
    for( int iFn   = 0; iFn   < nFineF  ; iFn++   ) {
    for( int iNX_F = 0; iNX_F < nDOFX_X1; iNX_F++ ) {
        k += 1;
        LX_X1[iNX_C][iFn][iNX_F] = LX_X1_Refined_Packed[k];
    }}}

    k = -1;
    for( int iNX_C = 0; iNX_C < nDOFX_X2; iNX_C++ ) {
    for( int iFn   = 0; iFn   < nFineF  ; iFn++   ) {
    for( int iNX_F = 0; iNX_F < nDOFX_X2; iNX_F++ ) {
        k += 1;
        LX_X2[iNX_C][iFn][iNX_F] = LX_X2_Refined_Packed[k];
    }}}

    k = -1;
    for( int iNX_C = 0; iNX_C < nDOFX_X3; iNX_C++ ) {
    for( int iFn   = 0; iFn   < nFineF  ; iFn++   ) {
    for( int iNX_F = 0; iNX_F < nDOFX_X3; iNX_F++ ) {
        k += 1;
        LX_X3[iNX_C][iFn][iNX_F] = LX_X3_Refined_Packed[k];
    }}}

    LX_X1_Up = LX_X1_Up_1D;
    LX_X1_Dn = LX_X1_Dn_1D;
    LX_X2_Up = LX_X2_Up_1D;
    LX_X2_Dn = LX_X2_Dn_1D;
    LX_X3_Up = LX_X3_Up_1D;
    LX_X3_Dn = LX_X3_Dn_1D;

    k = -1;
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
        k += 1;
        NodeNumberTableX[0][k] = iNX1;
        NodeNumberTableX[1][k] = iNX2;
        NodeNumberTableX[2][k] = iNX3;
        NodeNumberTableX3D[iNX1][iNX2][iNX3] = k;
    }}}

    /* iDimX == 0 */
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
        k = NodeNumberTableX3D[iNX1][iNX2][iNX3];
        NodeNumberTableX_X1[k] = k / nNodesX[0];
    }}}

    /* iDimX == 1 */
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      k = NodeNumberTableX3D[iNX1][iNX2][iNX3];
      NodeNumberTableX_X2[k] = k % nNodesX[0] + nNodesX[0] * iNX3;
    }}}

    /* iDimX == 2 */
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      k = NodeNumberTableX3D[iNX1][iNX2][iNX3];
      NodeNumberTableX_X3[k] = k % nDOFX_X1;
    }}}

    iGF_SqrtGm = iGF_SqtGm;
} /* END void InitializeMeshRefinement_DG */

void FinalizeMeshRefinement_DG()
{
    DeallocateArray( 3, NodeNumberTableX );
    DeallocateArray( nNodesX[0], nNodesX[1], NodeNumberTableX3D );
    DeallocateArray( NodeNumberTableX_X1 );
    DeallocateArray( NodeNumberTableX_X2 );
    DeallocateArray( NodeNumberTableX_X3 );
    DeallocateArray( LX_X3_Dn );
    DeallocateArray( LX_X3_Up );
    DeallocateArray( LX_X2_Dn );
    DeallocateArray( LX_X2_Dn );
    DeallocateArray( LX_X1_Up );
    DeallocateArray( LX_X1_Up );
    DeallocateArray( nDOFX_X3, nFineF, LX_X3 );
    DeallocateArray( nDOFX_X2, nFineF, LX_X2 );
    DeallocateArray( nDOFX_X1, nFineF, LX_X1 );
    DeallocateArray( nFineV, nDOFX, ProjectionMatrix_T );
    DeallocateArray( nFineV, nDOFX, ProjectionMatrix   );
    DeallocateArray( WeightsX_q );
    DeallocateArray( WeightsX_X3 );
    DeallocateArray( WeightsX_X2 );
    DeallocateArray( WeightsX_X1 );
} /* END void FinalizeMeshRefinement_DG */


/* PRIVATE FUNCTIONS */


void AllocateArray( int n0, Real * &A )
{
    A = new Real [n0];
}
void DeallocateArray( Real * &A )
{
    delete [] A; A = NULL;
}
void AllocateArray( int n0, int * &A )
{
    A = new int [n0];
}
void DeallocateArray( int * &A )
{
    delete [] A; A = NULL;
}

void AllocateArray( int n0, int n1, Real ** &A )
{
    A = new Real * [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new Real [n1];
    }
}
void DeallocateArray( int n0, Real ** &A )
{
    for( int i0 = 0; i0 < n0; i0++ ) {
        delete [] A[i0]; A[i0] = NULL;
    }
    delete [] A; A = NULL;
}
void AllocateArray( int n0, int n1, int ** &A )
{
    A = new int * [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new int [n1];
    }
}
void DeallocateArray( int n0, int ** &A )
{
    for( int i0 = 0; i0 < n0; i0++ ) {
        delete [] A[i0]; A[i0] = NULL;
    }
    delete [] A; A = NULL;
}

void AllocateArray( int n0, int n1, int n2, Real *** &A )
{
    A = new Real ** [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new Real * [n1];
        for( int i1 = 0; i1 < n1; i1++ ) {
            A[i0][i1] = new Real [n2];
        }
    }
}
void DeallocateArray( int n0, int n1, Real *** &A )
{
    for( int i0 = 0; i0 < n0; i0++ ) {
        for( int i1 = 0; i1 < n1; i1++ ) {
            delete [] A[i0][i1]; A[i0][i1] = NULL;
        }
        delete [] A[i0]; A[i0] = NULL;
    }
    delete [] A; A = NULL;
}
void AllocateArray( int n0, int n1, int n2, int *** &A )
{
    A = new int ** [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new int * [n1];
        for( int i1 = 0; i1 < n1; i1++ ) {
            A[i0][i1] = new int [n2];
        }
    }
}
void DeallocateArray( int n0, int n1, int *** &A )
{
    for( int i0 = 0; i0 < n0; i0++ ) {
        for( int i1 = 0; i1 < n1; i1++ ) {
            delete [] A[i0][i1]; A[i0][i1] = NULL;
        }
        delete [] A[i0]; A[i0] = NULL;
    }
    delete [] A; A = NULL;
}

} /* END namespace DG */
} /* END namespace amrex */
