#include <AMReX_FPhysBC.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

namespace
{
    // THIS MUST BE CONSISTENT WITH amrex_interpolater_module in AMReX_interpolater_mod.F90!!!
    static Array<Interpolater*> interp = {
        &amrex::pc_interp,               // 0 
        &amrex::node_bilinear_interp,    // 1
        &amrex::cell_bilinear_interp,    // 2
        &amrex::quadratic_interp,	 // 3
        &amrex::lincc_interp,	         // 4
        &amrex::cell_cons_interp,	 // 5
        &amrex::protected_interp,	 // 6
        &amrex::quartic_interp           // 7
    };
}

extern "C"
{
    void amrex_fi_fillpatch_single (MultiFab* mf, Real time, MultiFab* smf[], Real stime[], int ns,
				    int scomp, int dcomp, int ncomp, const Geometry* geom, 
                                    FPhysBC::fill_physbc_funptr_t fill)
    {
        FPhysBC pbc(fill);
	amrex::FillPatchSingleLevel(*mf, time, Array<MultiFab*>{smf, smf+ns}, 
				    Array<Real>{stime, stime+ns},
				    scomp, dcomp, ncomp, *geom, pbc);
    }

    void amrex_fi_fillpatch_two (MultiFab* mf, Real time,
				 MultiFab* cmf[], Real ct[], int nc,
				 MultiFab* fmf[], Real ft[], int nf,
				 int scomp, int dcomp, int ncomp,
				 const Geometry* cgeom, const Geometry* fgeom,
				 FPhysBC::fill_physbc_funptr_t cfill,
                                 FPhysBC::fill_physbc_funptr_t ffill,
				 int rr, int interp_id,
				 int* lo_bc[], int* hi_bc[])
    {
	Array<BCRec> bcs;
	// skip first scomp components
	for (int i = scomp; i < scomp+ncomp; ++i) {
	    bcs.emplace_back(lo_bc[i], hi_bc[i]);
	}

        FPhysBC cbc(cfill);
        FPhysBC fbc(ffill);
	amrex::FillPatchTwoLevels(*mf, time,
				  Array<MultiFab*>{cmf, cmf+nc}, Array<Real>{ct, ct+nc},
				  Array<MultiFab*>{fmf, fmf+nf}, Array<Real>{ft, ft+nf},
				  scomp, dcomp, ncomp,
				  *cgeom, *fgeom,
				  cbc, fbc,
				  IntVect{D_DECL(rr,rr,rr)},
				  interp[interp_id], bcs);
    }

    void amrex_fi_fillcoarsepatch (MultiFab* mf, Real time, const MultiFab* cmf,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry* cgeom, const Geometry* fgeom,
                                   FPhysBC::fill_physbc_funptr_t cfill,
                                   FPhysBC::fill_physbc_funptr_t ffill,
                                   int rr, int interp_id,
                                   int* lo_bc[], int* hi_bc[])
    {
	Array<BCRec> bcs;
        // skip first scomp components
	for (int i = scomp; i < scomp+ncomp; ++i) {
	    bcs.emplace_back(lo_bc[i], hi_bc[i]);
	}

        FPhysBC cbc(cfill);
        FPhysBC fbc(ffill);
        amrex::InterpFromCoarseLevel(*mf, time, *cmf, 
                                     scomp, dcomp, ncomp,
                                     *cgeom, *fgeom,
                                     cbc, fbc,
                                     IntVect{D_DECL(rr,rr,rr)},
                                     interp[interp_id], bcs);
    }
}
