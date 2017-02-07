#include <AMReX_MultiFabUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_average_down (const MultiFab* S_fine, MultiFab* S_crse,
			     const Geometry* fgeom, const Geometry* cgeom,
			     int scomp, int ncomp, int rr)
    {
	amrex::average_down(*S_fine, *S_crse, *fgeom, *cgeom, scomp, ncomp, rr);
    }
}
