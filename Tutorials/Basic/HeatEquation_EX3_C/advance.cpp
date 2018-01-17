
#include "myfunc.H"
#include "myfunc_F.H"

#include <AMReX_BCUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
              Real dt,
              const Geometry& geom,
              const BoxArray& grids, 
              const DistributionMapping& dmap, 
              const Vector<BCRec>& bc)
{

    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    // Fill non-periodic physical boundaries
    FillDomainBoundary(phi_old, geom, bc);

    // For MLMG solver
    int verbose = 2;
    int cg_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);

    const Real tol_rel  = 1.e-10;
    const Real tol_abs = 0.0;

    //
    // Implicit solve using MLABecLaplacian class
    // Assume composite solve, composite_solve = 1
    // 
    MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);
    
    mlabec.setMaxOrder(linop_maxorder);

    // Implement BCs 
    // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;

    for (int n = 0; n < phi_old.nComp(); ++n) 
    {
	for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
	    // lo-side BCs
	    if (bc[n].lo(idim) == BCType::int_dir) {
		bc_lo[idim] = LinOpBCType::Periodic;
	    }
	    else if (bc[n].lo(idim) == BCType::foextrap) {
		bc_lo[idim] = LinOpBCType::Neumann;
	    }
	    else if (bc[n].lo(idim) == BCType::ext_dir) {
		bc_lo[idim] = LinOpBCType::Dirichlet;
	    }
	    else {
		amrex::Abort("Invalid bc_lo");
	    }
	    
	    // hi-side BCs
	    if (bc[n].hi(idim) == BCType::int_dir) {
		bc_hi[idim] = LinOpBCType::Periodic;
	    }
	    else if (bc[n].hi(idim) == BCType::foextrap) {
		bc_hi[idim] = LinOpBCType::Neumann;
	    }
	    else if (bc[n].hi(idim) == BCType::ext_dir) {
		bc_hi[idim] = LinOpBCType::Dirichlet;
	    }
	    else {
		amrex::Abort("Invalid bc_hi");
	    }
	}
    }

    mlabec.setDomainBC(bc_lo, bc_hi);
    
    mlabec.setLevelBC(0, &phi_old);

    Real ascalar = 1.0;
    Real bscalar = 1.0;
    mlabec.setScalars(ascalar, bscalar);

    // Set up coefficient matrices
    MultiFab acoef(grids, dmap, 1, 0);
    MultiFab bcoef(grids, dmap, 1, 1);

    acoef.setVal(1.0);
    bcoef.setVal(dt);

    mlabec.setACoeffs(0, acoef);
	    
    std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
	const BoxArray& ba = amrex::convert(bcoef.boxArray(),
					    IntVect::TheDimensionVector(idim));
	face_bcoef[idim].define(ba, bcoef.DistributionMap(), 1, 0);

	face_bcoef[idim].setVal(dt);
    };

    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setCGVerbose(cg_verbose);

    // Solve linear system
    mlmg.solve({&phi_new}, {&phi_old}, tol_rel, tol_abs);
}

