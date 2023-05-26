#include <AMReX_DG.H>

#include <cstddef> /* For NULL */
#include <iostream> /* for std::cout/std::endl */

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>

namespace amrex::DG
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

    nFineV = (int)std::pow( 2, AMREX_SPACEDIM );
    nFineF = (int)std::pow( 2, AMREX_SPACEDIM-1 );

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

    AllocateArray( nNodesX[0], LX_X1_Up );
    AllocateArray( nNodesX[0], LX_X1_Dn );
    AllocateArray( nNodesX[1], LX_X2_Up );
    AllocateArray( nNodesX[1], LX_X2_Dn );
    AllocateArray( nNodesX[2], LX_X3_Up );
    AllocateArray( nNodesX[2], LX_X3_Dn );

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

    for( int iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
      LX_X1_Up[iNX1] = LX_X1_Up_1D[iNX1];
      LX_X1_Dn[iNX1] = LX_X1_Dn_1D[iNX1];
    }
    for( int iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
      LX_X2_Up[iNX2] = LX_X2_Up_1D[iNX2];
      LX_X2_Dn[iNX2] = LX_X2_Dn_1D[iNX2];
    }
    for( int iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
      LX_X3_Up[iNX3] = LX_X3_Up_1D[iNX3];
      LX_X3_Dn[iNX3] = LX_X3_Dn_1D[iNX3];
    }

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

#ifdef AMREX_USE_GPU

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &WeightsX_X1[0],
                      &WeightsX_X1[nDOFX_X1],
                      &WeightsX_X1[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &WeightsX_X2[0],
                      &WeightsX_X2[nDOFX_X2],
                      &WeightsX_X2[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &WeightsX_X3[0],
                      &WeightsX_X3[nDOFX_X3],
                      &WeightsX_X3[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &WeightsX_q[0],
                      &WeightsX_q[nDOFX],
                      &WeightsX_q[0]);

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &ProjectionMatrix[0],
                      &ProjectionMatrix[nFineV*nDOFX*nDOFX],
                      &ProjectionMatrix[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &ProjectionMatrix_T[0],
                      &ProjectionMatrix_T[nFineV*nDOFX*nDOFX],
                      &ProjectionMatrix_T[0]);

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X1[0],
                      &LX_X1[nDOFX_X1*nFineF*nDOFX_X1],
                      &LX_X1[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X2[0],
                      &LX_X2[nDOFX_X2*nFineF*nDOFX_X2],
                      &LX_X2[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X3[0],
                      &LX_X3[nDOFX_X3*nFineF*nDOFX_X3],
                      &LX_X3[0]);

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X1_Up[0],
                      &LX_X1_Up[nNodesX[0]],
                      &LX_X1_Up[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X1_Dn[0],
                      &LX_X1_Dn[nNodesX[0]],
                      &LX_X1_Dn[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X2_Up[0],
                      &LX_X2_Up[nNodesX[1]],
                      &LX_X2_Up[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X2_Dn[0],
                      &LX_X2_Dn[nNodesX[1]],
                      &LX_X2_Dn[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X3_Up[0],
                      &LX_X3_Up[nNodesX[2]],
                      &LX_X3_Up[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &LX_X3_Dn[0],
                      &LX_X3_Dn[nNodesX[2]],
                      &LX_X3_Dn[0]);

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &NodeNumberTableX[0],
                      &NodeNumberTableX[3*nDOFX],
                      &NodeNumberTableX[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &NodeNumberTableX3D[0],
                      &NodeNumberTableX3D[nDOFX],
                      &NodeNumberTableX3D[0]);

    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &NodeNumberTableX_X1[0],
                      &NodeNumberTableX_X1[nDOFX],
                      &NodeNumberTableX_X1[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &NodeNumberTableX_X2[0],
                      &NodeNumberTableX_X2[nDOFX],
                      &NodeNumberTableX_X2[0]);
    amrex::Gpu::copy( amrex::Gpu::hostToDevice,
                      &NodeNumberTableX_X3[0],
                      &NodeNumberTableX_X3[nDOFX],
                      &NodeNumberTableX_X3[0]);

#endif

} /* END void InitializeMeshRefinement_DG */

void FinalizeMeshRefinement_DG()
{
    DeallocateArray( NodeNumberTableX_X3 );
    DeallocateArray( NodeNumberTableX_X2 );
    DeallocateArray( NodeNumberTableX_X1 );
    DeallocateArray( nNodesX[0], nNodesX[1], NodeNumberTableX3D );
    DeallocateArray( 3, NodeNumberTableX );
    DeallocateArray( LX_X3_Dn );
    DeallocateArray( LX_X3_Up );
    DeallocateArray( LX_X2_Dn );
    DeallocateArray( LX_X2_Up );
    DeallocateArray( LX_X1_Dn );
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
#ifdef AMREX_USE_GPU
    A = (Real*)( The_Device_Arena() -> alloc( n0 * sizeof( Real ) ) );
#endif
}
void DeallocateArray( Real * &A )
{
#ifdef AMREX_USE_GPU
    The_Device_Arena() -> free( A );
#endif
    delete [] A; A = NULL;
}
void AllocateArray( int n0, int * &A )
{
    A = new int [n0];
#ifdef AMREX_USE_GPU
    A = (int*)( The_Device_Arena() -> alloc( n0 * sizeof( int ) ) );
#endif
}
void DeallocateArray( int * &A )
{
#ifdef AMREX_USE_GPU
    The_Device_Arena() -> free( A );
#endif
    delete [] A; A = NULL;
}

void AllocateArray( int n0, int n1, Real ** &A )
{
    A = new Real * [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new Real [n1];
    }
#ifdef AMREX_USE_GPU
    A = (Real**)( The_Device_Arena() -> alloc( n0*n1 * sizeof( Real ) ) );
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = (Real*)( The_Device_Arena() -> alloc( n1 * sizeof( Real ) ) );
    }
#endif
}
void DeallocateArray( int n0, Real ** &A )
{
#ifdef AMREX_USE_GPU
    for( int i0 = 0; i0 < n0; i0++ ) {
        The_Device_Arena() -> free( A[i0] );
    }
    The_Device_Arena() -> free( A );
#endif
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
#ifdef AMREX_USE_GPU
    A = (int**)( The_Device_Arena() -> alloc( n0*n1 * sizeof( int ) ) );
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = (int*)( The_Device_Arena() -> alloc( n1 * sizeof( int ) ) );
    }
#endif
}
void DeallocateArray( int n0, int ** &A )
{
#ifdef AMREX_USE_GPU
    for( int i0 = 0; i0 < n0; i0++ ) {
        The_Device_Arena() -> free( A[i0] );
    }
    The_Device_Arena() -> free( A );
#endif
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
#ifdef AMREX_USE_GPU
    A = (Real***)( The_Device_Arena() -> alloc( n0*n1*n2 * sizeof( Real ) ) );
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = (Real**)( The_Device_Arena()
                  -> alloc( n1*n2 * sizeof( Real ) ) );
        for( int i1 = 0; i1 < n1; i1++ ) {
            A[i0][i1] = (Real*)( The_Device_Arena()
                          -> alloc( n2 * sizeof( Real ) ) );
        }
    }
#endif
}
void DeallocateArray( int n0, int n1, Real *** &A )
{
#ifdef AMREX_USE_GPU
    for( int i0 = 0; i0 < n0; i0++ ) {
        for( int i1 = 0; i1 < n1; i1++ ) {
            The_Device_Arena() -> free( A[i0][i1] );
        }
        The_Device_Arena() -> free( A[i0] );
    }
    The_Device_Arena() -> free( A );
#endif
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
#ifdef AMREX_USE_GPU
    A = (int***)( The_Device_Arena() -> alloc( n0*n1*n2 * sizeof( int ) ) );
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = (int**)( The_Device_Arena()
                  -> alloc( n1*n2 * sizeof( int ) ) );
        for( int i1 = 0; i1 < n1; i1++ ) {
            A[i0][i1] = (int*)( The_Device_Arena()
                          -> alloc( n2 * sizeof( int ) ) );
        }
    }
#endif
}
void DeallocateArray( int n0, int n1, int *** &A )
{
#ifdef AMREX_USE_GPU
    for( int i0 = 0; i0 < n0; i0++ ) {
        for( int i1 = 0; i1 < n1; i1++ ) {
            The_Device_Arena() -> free( A[i0][i1] );
        }
        The_Device_Arena() -> free( A[i0] );
    }
    The_Device_Arena() -> free( A );
#endif
    for( int i0 = 0; i0 < n0; i0++ ) {
        for( int i1 = 0; i1 < n1; i1++ ) {
            delete [] A[i0][i1]; A[i0][i1] = NULL;
        }
        delete [] A[i0]; A[i0] = NULL;
    }
    delete [] A; A = NULL;
}

} /* END namespace amrex::DG */
