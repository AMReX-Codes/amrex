#include <AMReX_thornado.H>
#include <AMReX_REAL.H>

using namespace amrex;
using namespace thornado;

/* names of functions being defined must be all lowercase letters! */

extern "C"
{
  void amrex_fi_initializemeshrefinement_thornado
    ( int nNodes[], Real ProjectionMatrix[], Real WeightsX_q[] )
  {
    InitializeMeshRefinement_Thornado( nNodes, ProjectionMatrix, WeightsX_q );
  }

  void amrex_fi_finalizemeshrefinement_thornado()
  {
    FinalizeMeshRefinement_Thornado();
  }
}
