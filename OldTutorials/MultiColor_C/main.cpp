
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ForkJoin.H>

using namespace amrex;

namespace {
    int ncolors = 2;
    
    int ncomp  = 8;
    Real a     = 1.e-3;
    Real b     = 1.0;
    Real sigma = 10.0;  // controls the size of jump
    Real w     = 0.05;  // contols the width of the jump

    int  verbose = 0;
    Real tolerance_rel = 1.e-8;
    Real tolerance_abs = 0.0;
}

extern "C"
{
    void fort_set_rhs(double*, const int*, const int*, int,
		      const double*, double, double, double, double);
    void fort_set_coef(double*, const int*, const int*, 
		       double*, const int*, const int*, 
		       int, const double*, double, double);
}

void setup_rhs(MultiFab& rhs, const Geometry& geom);
void setup_coeffs(MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom);
void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const Vector<MultiFab*>& beta,
	   const Geometry& geom);
void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const Vector<MultiFab*>& beta, 
		   const Geometry& geom);
void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const Vector<MultiFab*>& beta, 
			    const Geometry& geom);

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BoxArray ba;
    Geometry geom;
    {
	ParmParse pp;
	
        pp.query("ncolors", ncolors);

	pp.query("verbose", verbose);

	int n_cell, max_grid_size;
	pp.get("n_cell", n_cell);
        if (pp.contains("ncomp")) {
            pp.get("ncomp", ncomp);
        }
	pp.get("max_grid_size", max_grid_size);

	Box domain(IntVect(AMREX_D_DECL(       0,       0,       0)),
		   IntVect(AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)));

	ba.define(domain);
	ba.maxSize(max_grid_size);

	RealBox real_box;
	for (int n = 0; n < BL_SPACEDIM; n++) {
	    real_box.setLo(n, 0.0);
	    real_box.setHi(n, 1.0);
	}

	int coord = 0;
 
	geom.define(domain,&real_box,coord);
    }

    DistributionMapping dm{ba};

    MultiFab rhs(ba, dm, ncomp, 0);
    setup_rhs(rhs, geom);

    MultiFab alpha(ba, dm, ncomp, 0);
    Vector<MultiFab> beta(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        beta[i].define(amrex::convert(ba, IntVect::TheDimensionVector(i)),
                       dm, ncomp, 0);
    }
    setup_coeffs(alpha, amrex::GetVecOfPtrs(beta), geom);

    MultiFab soln(ba, dm, ncomp, 0);

    solve(soln, rhs, alpha, amrex::GetVecOfPtrs(beta), geom);

    VisMF::Write(soln, "soln");

    amrex::Finalize();
}

void setup_rhs(MultiFab& rhs, const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
    {
	const int* rlo = rhs[mfi].loVect();
	const int* rhi = rhs[mfi].hiVect();

	fort_set_rhs(rhs[mfi].dataPtr(),rlo, rhi, rhs.nComp(),
		     dx, a, b, sigma, w);
    }
}

void setup_coeffs(MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    alpha.setVal(1.0);

#if (BL_SPACEDIM == 3)
    amrex::Abort("2D only");
#endif

    for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) 
    {
	FArrayBox& betax = (*beta[0])[mfi];
	FArrayBox& betay = (*beta[1])[mfi];

	fort_set_coef(betax.dataPtr(), betax.loVect(), betax.hiVect(),
		      betay.dataPtr(), betay.loVect(), betay.hiVect(),
		      beta[0]->nComp(), dx, sigma, w);
    }
}

void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom)
{
    // evenly split ranks among ncolors tasks
    ForkJoin fj(ncolors);

    // register how to copy multifabs to/from tasks
    fj.reg_mf    (rhs  , "rhs"  , ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf    (alpha, "alpha", ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf_vec(beta , "beta" , ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf    (soln , "soln" , ForkJoin::Strategy::split, ForkJoin::Intent::out);

    Vector<ForkJoin::ComponentSet> comp_split {{0,5}, {5,8}, {8,12}, {12,16}};

    fj.modify_split("rhs", 0, comp_split);
    fj.modify_split("alpha", 0, comp_split);
    for (int i = 0; i < beta.size(); ++i) {
        fj.modify_split("beta", i, comp_split);
    }
    fj.modify_split("soln", 0, comp_split);

    // can reuse ForkJoin object for multiple fork-join invocations
    // creates forked multifabs only first time around, reuses them thereafter
    for (int i = 0; i < 2; ++i) {
        // issue fork-join
        fj.fork_join(
            [&geom] (ForkJoin &f) {
                colored_solve(f.get_mf("soln"), f.get_mf("rhs"), f.get_mf("alpha"),
                              f.get_mf_vec("beta"), geom);
            }
        );
    }
}

void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const Vector<MultiFab*>& beta, 
		   const Geometry& geom)
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();

    if (rhs.nComp() == 1)
    {
	single_component_solve(soln, rhs, alpha, beta, geom);
    }
    else
    {	    
	for (int i = 0; i < soln.nComp(); ++i) {
	    MultiFab ssoln (ba, dm, 1, 1);
            ssoln.setVal(0.0);

	    MultiFab srhs  (rhs, amrex::make_alias, i, 1);
	    MultiFab salpha(alpha, amrex::make_alias, i, 1);
            Vector<std::unique_ptr<MultiFab> > sbeta(AMREX_SPACEDIM);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                sbeta[idim].reset(new MultiFab(*beta[idim], amrex::make_alias, i, 1));
            }

	    single_component_solve(ssoln, srhs, salpha,
				   amrex::GetVecOfPtrs(sbeta), geom);
	    
	    MultiFab::Copy(soln, ssoln, 0, i, 1, 0);
	}
    }
}

void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const Vector<MultiFab*>& beta, 
			    const Geometry& geom)
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();

    MLABecLaplacian mlabec({geom}, {ba}, {dm});

    mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)},
                       {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)});
    mlabec.setLevelBC(0, &soln);

    mlabec.setScalars(a, b);
    mlabec.setACoeffs(0, alpha);
    mlabec.setBCoeffs(0, {AMREX_D_DECL(beta[0], beta[1], beta[2])});

    MLMG mlmg(mlabec);
    mlmg.setVerbose(verbose);
    mlmg.solve({&soln}, {&rhs}, tolerance_rel, tolerance_abs);
}

