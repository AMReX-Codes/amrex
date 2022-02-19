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
extern int nFine;
extern Real *** ProjectionMatrix;
extern Real *** ProjectionMatrix_T;
extern Real * WeightsX_q;
extern Real VolumeRatio;

void InitializeMeshRefinement_Thornado
       ( int N[], Real ProjMatrix[], Real WeightsX[] );
void FinalizeMeshRefinement_Thornado();

void AllocateArray( int n0, Real * &A );
void DeallocateArray( Real * &A );

void AllocateArray( int n0, int n1, int n2, Real *** &A );
void DeallocateArray( int n0, int n1, Real *** &A );

} /* END namespace thornado */
} /* END namespace amrex */

#endif // AMREX_THORNADO_H_
