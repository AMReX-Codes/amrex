#include <AMReX_AlgMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLALaplacian.H>
#include <AMReX_MLPoisson.H>

namespace amrex {

// Single precision copies used by mixed precision, in their own file so
// that they compile in parallel with AMReX_MLMG.cpp.
#ifdef AMREX_MLMG_HAS_FLOAT
template class MLMGT<fMultiFab>;
template class MLLinOpT<fMultiFab>;
template class MLCellLinOpT<fMultiFab>;
template class MLCellABecLapT<fMultiFab>;
template class MLABecLaplacianT<fMultiFab>;
template class MLALaplacianT<fMultiFab>;
template class MLPoissonT<fMultiFab>;
#endif

}
