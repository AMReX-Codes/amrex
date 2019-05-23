#include <AMReX_FPhysBC.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

namespace
{
    // THIS MUST BE CONSISTENT WITH amrex_interpolater_module in AMReX_interpolater_mod.F90!!!
    static Vector<Interpolater*> interp = {
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

namespace {
    
    typedef void (*INTERP_HOOK) (const int* lo, const int*hi,
                                 Real* d, const int* dlo, const int* dhi, const int nd,
                                 const int icomp, const int ncomp);

    class FIInterpHook final
        : public InterpHook
    {
    public:
        explicit FIInterpHook (INTERP_HOOK a_f) : m_f(a_f) {}
        virtual void operator() (FArrayBox& fab, const Box& bx, int icomp, int ncomp) const final
        {
            if (m_f) {
                m_f(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(fab), fab.nComp(),
                    icomp+1, ncomp);
                // m_f is a fortran function expecting 1-based index
            }
        }
    private:
        INTERP_HOOK m_f;
    };
}

extern "C"
{
    void amrex_fi_fillpatch_single (MultiFab* mf, Real time, MultiFab* smf[], Real stime[], int ns,
				    int scomp, int dcomp, int ncomp, const Geometry* geom, 
                                    FPhysBC::fill_physbc_funptr_t fill)
    {
        FPhysBC pbc(fill, geom);
	amrex::FillPatchSingleLevel(*mf, time, Vector<MultiFab*>{smf, smf+ns}, 
				    Vector<Real>{stime, stime+ns},
				    scomp, dcomp, ncomp, *geom, pbc, 0);
    }

    void amrex_fi_fillpatch_two (MultiFab* mf, Real time,
				 MultiFab* cmf[], Real ct[], int nc,
				 MultiFab* fmf[], Real ft[], int nf,
				 int scomp, int dcomp, int ncomp,
				 const Geometry* cgeom, const Geometry* fgeom,
				 FPhysBC::fill_physbc_funptr_t cfill,
                                 FPhysBC::fill_physbc_funptr_t ffill,
				 int rr, int interp_id,
				 int* lo_bc[], int* hi_bc[],
                                 INTERP_HOOK pre_interp, INTERP_HOOK post_interp)
    {
	Vector<BCRec> bcs(ncomp);
	for (int i = 0; i < ncomp; ++i) {
	    bcs.emplace_back(lo_bc[i], hi_bc[i]);
	}

        FPhysBC cbc(cfill, cgeom);
        FPhysBC fbc(ffill, fgeom);
	amrex::FillPatchTwoLevels(*mf, time,
				  Vector<MultiFab*>{cmf, cmf+nc}, Vector<Real>{ct, ct+nc},
				  Vector<MultiFab*>{fmf, fmf+nf}, Vector<Real>{ft, ft+nf},
				  scomp, dcomp, ncomp,
				  *cgeom, *fgeom,
				  cbc, 0, fbc, 0,
				  IntVect{AMREX_D_DECL(rr,rr,rr)},
				  interp[interp_id], bcs, 0,
                                  FIInterpHook(pre_interp),
                                  FIInterpHook(post_interp));
    }

    void amrex_fi_fillcoarsepatch (MultiFab* mf, Real time, const MultiFab* cmf,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry* cgeom, const Geometry* fgeom,
                                   FPhysBC::fill_physbc_funptr_t cfill,
                                   FPhysBC::fill_physbc_funptr_t ffill,
                                   int rr, int interp_id,
                                   int* lo_bc[], int* hi_bc[],
                                   INTERP_HOOK pre_interp, INTERP_HOOK post_interp)
    {
	Vector<BCRec> bcs(ncomp);
	for (int i = 0; i < ncomp; ++i) {
	    bcs.emplace_back(lo_bc[i], hi_bc[i]);
	}

        FPhysBC cbc(cfill, cgeom);
        FPhysBC fbc(ffill, fgeom);
        amrex::InterpFromCoarseLevel(*mf, time, *cmf, 
                                     scomp, dcomp, ncomp,
                                     *cgeom, *fgeom,
                                     cbc, 0, fbc, 0,
                                     IntVect{AMREX_D_DECL(rr,rr,rr)},
                                     interp[interp_id], bcs, 0,
                                     FIInterpHook(pre_interp),
                                     FIInterpHook(post_interp));
    }
}
