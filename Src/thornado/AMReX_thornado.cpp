#include <AMReX_thornado.H>

#include <cstddef> /* For NULL */
#include <iostream> /* for std::cout/std::endl */
#include <cmath> /* for pow */

#include <AMReX_REAL.H>

namespace amrex
{
namespace thornado
{

int nNodes[3];
int nDOFX;
int nFine;
Real *** ProjectionMatrix   = NULL;
Real *** ProjectionMatrix_T = NULL;
Real * WeightsX_q = NULL;

void InitializeMeshRefinement_Thornado
       ( int N[], Real ProjMatrix[], Real WeightsX[] )
{
    nNodes[0] = N[0];
    nNodes[1] = N[1];
    nNodes[2] = N[2];

    nFine = pow( 2, AMREX_SPACEDIM );

    nDOFX = N[0] * N[1] * N[2];

    AllocateArray( nFine, nDOFX, nDOFX, ProjectionMatrix   );
    AllocateArray( nFine, nDOFX, nDOFX, ProjectionMatrix_T );
    AllocateArray( nDOFX, WeightsX_q );

    int k = -1;
    for( int iFine = 0; iFine < nFine; iFine++ ) {
    for( int jNX   = 0; jNX   < nDOFX; jNX++   ) {
    for( int iNX   = 0; iNX   < nDOFX; iNX++   ) {
        k += 1;
        ProjectionMatrix  [iFine][iNX][jNX] = ProjMatrix[k];
        ProjectionMatrix_T[iFine][jNX][iNX] = ProjMatrix[k];
    }}}

    for( int iNX = 0; iNX < nDOFX; iNX++ ) {
      WeightsX_q[iNX] = WeightsX[iNX];
    }

} /* END void InitializeMeshRefinement_Thornado */

void FinalizeMeshRefinement_Thornado()
{
    DeallocateArray( WeightsX_q );
    DeallocateArray( nFine, nDOFX, ProjectionMatrix_T );
    DeallocateArray( nFine, nDOFX, ProjectionMatrix   );

} /* END void FinalizeMeshRefinement_Thornado */


/* PRIVATE FUNCTIONS */


void AllocateArray( int n0, Real * &A )
{
    A = new Real [n0];
} /* END void AallocateArray */

void DeallocateArray( Real * &A )
{
    delete [] A; A = NULL;
} /* END void DeallocateArray */

void AllocateArray( int n0, int n1, int n2, Real *** &A )
{
    A = new Real ** [n0];
    for( int i0 = 0; i0 < n0; i0++ ) {
        A[i0] = new Real * [n1];
        for( int i1 = 0; i1 < n1; i1++ ) {
            A[i0][i1] = new Real [n2];
        }
    }
} /* END void AallocateArray */

void DeallocateArray( int n0, int n1, Real *** &A )
{
    for( int i0 = 0; i0 < n0; i0++ ) {
    for( int i1 = 0; i1 < n1; i1++ ) {
        delete [] A[i0][i1]; A[i0][i1] = NULL; }
        delete [] A[i0]; A[i0] = NULL; }
    delete [] A; A = NULL;
} /* END void DeallocateArray */

} /* END namespace thornado */
} /* END namespace amrex */
