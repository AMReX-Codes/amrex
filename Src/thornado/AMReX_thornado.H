#ifndef AMREX_THORNADO_H_
#define AMREX_THORNADO_H_

#include <AMReX_REAL.H>

namespace amrex
{
namespace thornado
{

const Real Zero = 0.0;
const Real One  = 1.0;

extern int nNodes[3];
extern int nDOFX;
extern int nDOFX_X1;
extern int nDOFX_X2;
extern int nDOFX_X3;
extern int nFineV;
extern int nFineF;
extern Real *** ProjectionMatrix;
extern Real *** ProjectionMatrix_T;
extern Real *** LX_X1;
extern Real *** LX_X2;
extern Real *** LX_X3;
extern Real * WeightsX_q;
extern Real * WeightsX_X1;
extern Real * WeightsX_X2;
extern Real * WeightsX_X3;
extern Real VolumeRatio;
extern Real FaceRatio;

void InitializeMeshRefinement_Thornado
       ( int N[], Real ProjMatrix[],
         Real WeightsX1[], Real WeightsX2[], Real WeightsX3[],
         Real LX_X1_Refined_Packed[],
         Real LX_X2_Refined_Packed[],
         Real LX_X3_Refined_Packed[] );

void FinalizeMeshRefinement_Thornado();

void AllocateArray( int n0, Real * &A );
void DeallocateArray( Real * &A );

void AllocateArray( int n0, int n1, int n2, Real *** &A );
void DeallocateArray( int n0, int n1, Real *** &A );

} /* END namespace thornado */
} /* END namespace amrex */

#endif // AMREX_THORNADO_H_
