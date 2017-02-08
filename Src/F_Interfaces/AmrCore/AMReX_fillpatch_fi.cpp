#include <AMReX_FPhysBC.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_fillpatch_single (MultiFab* mf, Real time, MultiFab* smf[], Real stime[], 
				    int scomp, int dcomp, int ncomp, 
				    const Geometry* geom, FPhysBC* pbc, int n)
    {
	amrex::FillPatchSingleLevel(*mf, time, Array<MultiFab*>{smf, smf+n}, 
				    Array<Real>{stime, stime+n},
				    scomp, dcomp, ncomp, *geom, *pbc);
    }
}
