#include <AMReX_thornado.H>
#include <AMReX_REAL.H>

using namespace amrex;
using namespace thornado;

/* names of declared functions must be all lowercase! */

extern "C"
{
  void amrex_fi_initializemeshrefinement_thornado
    ( int nNodes[], Real ProjectionMatrix[], Real WeightsX_q[] )
  {
    InitializeMeshRefinement_Thornado( nNodes, ProjectionMatrix, WeightsX_q );
  }
}
