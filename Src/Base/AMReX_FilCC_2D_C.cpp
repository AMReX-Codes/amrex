#include <AMReX_FilCC_C.H>

namespace amrex {

AMREX_GPU_DEVICE
void
filcc_cell (const IntVect& iv, FArrayBox& dest,
            const int dcomp, const int numcomp,
            GeometryData const& geom, const Real time,
            const BCRec* bcr, const int bcomp,
            const int orig_comp)
{

}

}
