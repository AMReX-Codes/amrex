#ifndef AMREX_DG_H_
#define AMREX_DG_H_

#include <AMReX_REAL.H>

namespace amrex::DG
{

const Real Zero = 0.0;
const Real One  = 1.0;

extern int nFineV;
extern int nFineF;

extern Real VolumeRatio;
extern Real FaceRatio;

extern int nNodes[3];
extern int nDOFX;
extern int nDOFX_X1;
extern int nDOFX_X2;
extern int nDOFX_X3;

extern Real * WeightsX_X1;
extern Real * WeightsX_X2;
extern Real * WeightsX_X3;
extern Real * WeightsX_q;

extern Real *** ProjectionMatrix;
extern Real *** ProjectionMatrix_T;

extern Real *** LX_X1;
extern Real *** LX_X2;
extern Real *** LX_X3;

extern Real * LX_X1_Up;
extern Real * LX_X1_Dn;
extern Real * LX_X2_Up;
extern Real * LX_X2_Dn;
extern Real * LX_X3_Up;
extern Real * LX_X3_Dn;

extern int **  NodeNumberTableX;
extern int *** NodeNumberTableX3D;
extern int *   NodeNumberTableX_X1;
extern int *   NodeNumberTableX_X2;
extern int *   NodeNumberTableX_X3;

extern int iGF_SqrtGm;

void InitializeMeshRefinement_DG
       ( const int N[],
         const Real ProjMatrix[],
         const Real WeightsX1[],
         const Real WeightsX2[],
         const Real WeightsX3[],
         const Real LX_X1_Refined_Packed[],
         const Real LX_X2_Refined_Packed[],
         const Real LX_X3_Refined_Packed[],
         const Real LX_X1_Up_1D[],
         const Real LX_X1_Dn_1D[],
         const Real LX_X2_Up_1D[],
         const Real LX_X2_Dn_1D[],
         const Real LX_X3_Up_1D[],
         const Real LX_X3_Dn_1D[],
         int iGF_SqtGm );

void FinalizeMeshRefinement_DG();

void AllocateArray  ( int n0, Real * &A );
void DeallocateArray(         Real * &A );
void AllocateArray  ( int n0, int * &A );
void DeallocateArray(         int * &A );

void AllocateArray  ( int n0, int n1, Real ** &A );
void DeallocateArray( int n0,         Real ** &A );
void AllocateArray  ( int n0, int n1, int ** &A );
void DeallocateArray( int n0,         int ** &A );

void AllocateArray  ( int n0, int n1, int n2, Real *** &A );
void DeallocateArray( int n0, int n1,         Real *** &A );
void AllocateArray  ( int n0, int n1, int n2, int *** &A );
void DeallocateArray( int n0, int n1,         int *** &A );

} /* END namespace amrex::DG */

#endif // AMREX_DG_H_
