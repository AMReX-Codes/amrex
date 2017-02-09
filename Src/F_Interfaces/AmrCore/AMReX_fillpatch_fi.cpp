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

    void amrex_fi_fillpatch_two (int interp_id)
    {
	// THIS MUST BE CONSISTENT WITH amrex_interpolater_module in AMReX_interpolater_mod.F90!!!
	static Array<Interpolater*> interp = {
	    &amrex::pc_interp,               // 0 
	    &amrex::node_bilinear_interp,    // 1
	    &amrex::cell_bilinear_interp,    // 2
	    &amrex::quadratic_interp,	     // 3
	    &amrex::lincc_interp,	     // 4
	    &amrex::cell_cons_interp,	     // 5
	    &amrex::protected_interp,	     // 6
	    &amrex::quartic_interp           // 7
	};
    }
}
