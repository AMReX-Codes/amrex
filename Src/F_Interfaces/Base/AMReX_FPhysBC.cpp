#include <AMReX_FPhysBC.H>

using namespace amrex;

void
amrex::FPhysBC::FillBoundary (MultiFab& mf, int scomp, int ncomp, Real time)
{
    if (fill_physbc != nullptr) {
	fill_physbc(&mf, scomp+1, ncomp, time, geom);
    } else {
	amrex::Abort("FPhysBC::fill_physbc is null");
    }
}
