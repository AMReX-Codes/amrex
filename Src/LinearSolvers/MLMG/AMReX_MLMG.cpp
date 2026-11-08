#include <AMReX_AlgMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLALaplacian.H>
#include <AMReX_MLPoisson.H>

namespace amrex {

// Explicit instantiation keeps applications from compiling these templates
// and works around nvcc issues.  The single precision versions are in
// AMReX_MLMG_float.cpp.

template class MLMGT<MultiFab>;
template class MLLinOpT<MultiFab>;
template class MLCellLinOpT<MultiFab>;
template class MLCellABecLapT<MultiFab>;
template class MLABecLaplacianT<MultiFab>;
template class MLALaplacianT<MultiFab>;
template class MLPoissonT<MultiFab>;

}
