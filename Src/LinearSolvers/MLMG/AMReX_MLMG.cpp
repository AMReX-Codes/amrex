#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLALaplacian.H>
#include <AMReX_MLPoisson.H>

namespace amrex {

// Explicit instantiation of these templates are necessary for working
// around nvcc issues.

template class MLMGT<MultiFab>;
template class MLLinOpT<MultiFab>;
template class MLCellLinOpT<MultiFab>;
template class MLCellABecLapT<MultiFab>;
template class MLABecLaplacianT<MultiFab>;
template class MLALaplacianT<MultiFab>;
template class MLPoissonT<MultiFab>;

}
