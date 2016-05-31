
#include <MultiFab.H>
#include <ParmParse.H>
#include <VisMF.H>
#include <Geometry.H>
#include <BndryData.H>
#include <LO_BCTYPES.H>
#include <MultiGrid.H>

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
void setup_coeffs(MultiFab& alpha, PArray<MultiFab>& beta, const Geometry& geom);
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);
void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const PArray<MultiFab>& beta,
	   const Geometry& geom);
void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const PArray<MultiFab>& beta, 
		   const Geometry& geom);
void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const PArray<MultiFab>& beta, 
			    const Geometry& geom);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BoxArray ba;
    Geometry geom;
    {
	ParmParse pp;
	
	pp.query("verbose", verbose);

	int n_cell, max_grid_size;
	pp.get("n_cell", n_cell);
	pp.get("max_grid_size", max_grid_size);

	Box domain(IntVect(D_DECL(       0,       0,       0)),
		   IntVect(D_DECL(n_cell-1,n_cell-1,n_cell-1)));

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

    MultiFab rhs(ba, ncomp, 0);
    setup_rhs(rhs, geom);

    MultiFab alpha(ba, ncomp, 0);
    PArray<MultiFab> beta(BL_SPACEDIM, PArrayManage);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	beta.set(i, new MultiFab(ba, ncomp, 0, Fab_allocate, IntVect::TheDimensionVector(i)));
    }
    setup_coeffs(alpha, beta, geom);

    MultiFab soln(ba, ncomp, 0);

    solve(soln, rhs, alpha, beta, geom);

    VisMF::Write(soln, "soln");

    BoxLib::Finalize();
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

void setup_coeffs(MultiFab& alpha, PArray<MultiFab>& beta, const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    alpha.setVal(1.0);

#if (BL_SPACEDIM == 3)
    BoxLib::Abort("2D only");
#endif

    for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) 
    {
	FArrayBox& betax = beta[0][mfi];
	FArrayBox& betay = beta[1][mfi];

	fort_set_coef(betax.dataPtr(), betax.loVect(), betax.hiVect(),
		      betay.dataPtr(), betay.loVect(), betay.hiVect(),
		      beta[0].nComp(), dx, sigma, w);
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
	   const MultiFab& alpha, const PArray<MultiFab>& beta, const Geometry& geom)
{
    int ncolors = ParallelDescriptor::NColors();
    int cncomp = ncomp / ncolors;
    BL_ASSERT(cncomp*ncolors == ncomp);

    const BoxArray& ba = rhs.boxArray();

    PArray<MultiFab> csoln(ncolors, PArrayManage);
    PArray<MultiFab> crhs(ncolors, PArrayManage);
    PArray<MultiFab> calpha(ncolors, PArrayManage);
    Array<PArray<MultiFab> > cbeta(ncolors);

    for (int i = 0; i < ncolors; ++i) {
	ParallelDescriptor::Color color = ParallelDescriptor::Color(i);

	// Build colored MFs.
	csoln.set(i, new MultiFab(ba, cncomp, 0, color));
	crhs.set(i, new MultiFab(ba, cncomp, 0, color));
	calpha.set(i, new MultiFab(ba, cncomp, 0, color));
	cbeta[i].resize(BL_SPACEDIM, PArrayManage);
	for (int j = 0; j < BL_SPACEDIM; ++j) {
	    cbeta[i].set(j, new MultiFab(beta[j].boxArray(), cncomp, 0, color));
	}

	// Parallel copy data into colored MFs.
	crhs[i].copy(rhs, cncomp*i, 0, cncomp);
	calpha[i].copy(alpha, cncomp*i, 0, cncomp);
	for (int j = 0; j < BL_SPACEDIM; ++j) {
	    cbeta[i][j].copy(beta[j], cncomp*i, 0, cncomp);
	}
    }

    int icolor = ParallelDescriptor::SubCommColor().to_int();
    colored_solve(csoln[icolor], crhs[icolor], calpha[icolor], cbeta[icolor], geom);

    // Copy solution back from colored MFs.
    for (int i = 0; i < ncolors; ++i) {
	soln.copy(csoln[i], 0, cncomp*i, cncomp);
    }
}

void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const PArray<MultiFab>& beta, 
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
	    MultiFab ssoln(ba, 1, 1, dm);
	    MultiFab srhs(ba, 1, 0, dm);
	    MultiFab salpha(ba, 1, 0, dm);
	    PArray<MultiFab> sbeta(BL_SPACEDIM, PArrayManage);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		sbeta.set(j, new MultiFab(beta[j].boxArray(), 1, 0, dm));
	    }
	    
	    MultiFab::Copy(ssoln , soln , i, 0, 1, 0);
	    MultiFab::Copy(srhs  , rhs  , i, 0, 1, 0);
	    MultiFab::Copy(salpha, alpha, i, 0, 1, 0);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		MultiFab::Copy(sbeta[j], beta[j], i, 0, 1, 0);
	    }
	    
	    single_component_solve(ssoln, srhs, salpha, sbeta, geom);
	    
	    MultiFab::Copy(soln, ssoln, 0, i, 1, 0);
	}
    }
}

void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const PArray<MultiFab>& beta, 
			    const Geometry& geom)
{
    const ParallelDescriptor::Color color = soln.color();
    const BoxArray& ba = soln.boxArray();
    const Real* dx = geom.CellSize();

    BndryData bd(ba, 1, geom, color);
    set_boundary(bd, rhs, 0);

    ABecLaplacian abec_operator(bd, dx);
    abec_operator.setScalars(a, b);
    abec_operator.setCoefficients(alpha, beta);

    MultiGrid mg(abec_operator);
    mg.setVerbose(verbose);

    soln.setVal(0.0);

    mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
}

