
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_Geometry.H>
#include <AMReX_BndryData.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ForkJoin.H>

using namespace amrex;

namespace {
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
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);
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
	
	pp.query("verbose", verbose);

	int n_cell, max_grid_size;
	pp.get("n_cell", n_cell);
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
    Vector<std::unique_ptr<MultiFab> > beta(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        beta[i].reset(new MultiFab(amrex::convert(ba, IntVect::TheDimensionVector(i)),
                                   dm, ncomp, 0));
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

// Dirichlet only in this test
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp)
{
    Real bc_value = 0.0;

    const Geometry& geom = bd.getGeom();
    const Real* dx = geom.CellSize();

    for (int n=0; n<BL_SPACEDIM; ++n) {
	for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
	    int i = mfi.index(); 
	    
	    const Box& bx = mfi.validbox();
	    
	    // Our default will be that the face of this grid is either touching another grid
	    //  across an interior boundary or a periodic boundary.  We will test for the other
	    //  cases below.
	    {
		// Define the type of boundary conditions to be Dirichlet (even for periodic)
		bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
		bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	
		// Set the boundary conditions to the cell centers outside the domain
		bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
		bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
	    }
	    
	    // Now test to see if we should override the above with Dirichlet physical bc's

	    // We are on the low side of the domain in coordinate direction n
	    if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
		// Set the boundary conditions to live exactly on the faces of the domain
		bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  
		// Set the Dirichlet/Neumann boundary values 
		bd.setValue(Orientation(n, Orientation::low) ,i, bc_value);
	  
		// Define the type of boundary conditions 
		bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	    }
	
	    // We are on the high side of the domain in coordinate direction n
	    if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
		// Set the boundary conditions to live exactly on the faces of the domain
		bd.setBoundLoc(Orientation(n, Orientation::high) ,i,0.0 );
		
		// Set the Dirichlet/Neumann boundary values
		bd.setValue(Orientation(n, Orientation::high) ,i, bc_value);
		
		// Define the type of boundary conditions 
		bd.setBoundCond(Orientation(n, Orientation::high) ,i,comp,LO_DIRICHLET);
	    }
	}
    }
}

void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom)
{
    // evenly split ranks among NColors() tasks
    ForkJoin fj(ParallelDescriptor::NColors());

    // register how to copy multifabs to/from tasks
    fj.reg_mf    (rhs  , "rhs"  , ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf    (alpha, "alpha", ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf_vec(beta , "beta" , ForkJoin::Strategy::split, ForkJoin::Intent::in);
    fj.reg_mf    (soln , "soln" , ForkJoin::Strategy::split, ForkJoin::Intent::out);

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
	    MultiFab srhs  (ba, dm, 1, 0);
	    MultiFab salpha(ba, dm, 1, 0);
	    Vector<std::unique_ptr<MultiFab> > sbeta(BL_SPACEDIM);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		sbeta[j].reset(new MultiFab(beta[j]->boxArray(), dm, 1, 0));
	    }
	    
	    MultiFab::Copy(ssoln , soln , i, 0, 1, 0);
	    MultiFab::Copy(srhs  , rhs  , i, 0, 1, 0);
	    MultiFab::Copy(salpha, alpha, i, 0, 1, 0);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		MultiFab::Copy(*sbeta[j], *beta[j], i, 0, 1, 0);
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
    const Real* dx = geom.CellSize();

    BndryData bd(ba, dm, 1, geom);
    set_boundary(bd, rhs, 0);

    ABecLaplacian abec_operator(bd, dx);
    abec_operator.setScalars(a, b);
    abec_operator.setCoefficients(alpha, beta);

    MultiGrid mg(abec_operator);
    mg.setVerbose(verbose);

    soln.setVal(0.0);

    mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
}

