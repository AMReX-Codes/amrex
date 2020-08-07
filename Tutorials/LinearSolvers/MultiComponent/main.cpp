//
// Tutorial:    MultiComponent Linear Solve
//
// File:        main.cpp
//
// Author:      Brandon Runnels
//              University of Colorado Colorado Springs
//              brunnels@uccs.edu
//              solids.uccs.edu
//
// Date:        September 3, 2019
// 
// Description: This tutorial demonstrates how to implement a multi-component
//              nodal linear operator. This tutorial demonstrates the
//              "CFStrategy::ghostnodes" method for computing the reflux at
//              the coarse/fine boundary. 
//
// See also:    ./MCNodalLinOp.H, ./MCNodalLinOp.cpp
//              for implementation of the MC linear operator.
// 

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Machine.H>
#include <AMReX_MLMG.H>
#include <AMReX_PlotFileUtil.H>
#include "MCNodalLinOp.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    Initialize(argc, argv);
    {
    
    //
    // Read in mesh structure (# amr levels, # nodes)
    // Geometry will be a nnodes x nnodes (x nnodes) grid
    // refined in the center up to nlevels.
    //
    struct {
        int nlevels = 3;
        int nnodes = 32;
        int max_grid_size = 10000000;
    } mesh;
    {
        ParmParse pp("mesh");
        pp.query("nlevels",mesh.nlevels);
        pp.query("nnodes",mesh.nnodes);
        pp.query("max_grid_size",mesh.max_grid_size);
    }
    
    //
    // Read in linear operator parameters:
    //    ncomp = number of components for a MC solve
    //    coeff = ncomp x ncomp list of coefficients
    //
    struct {
        int ncomp=1;
        Vector<Real> coeff = {1.0};
    } op;
    {
        ParmParse pp("op");
        pp.query("ncomp",op.ncomp);
        pp.queryarr("coeff",op.coeff);
    }

    //
    // Read in MLMG solver parameters
    //
    struct {
        int verbose = -1;
        int cg_verbose = -1;
        int max_iter = -1;
        int fixed_iter = -1;
        int max_fmg_iter = -1;
        int linop_maxorder = -1;
        int agglomeration = -1;
        int consolidation = -1;
        int max_coarsening_level = -1 ;
    } mlmg;
    {
        ParmParse pp("mlmg");
        pp.query("verbose",mlmg.verbose);
        pp.query("cg_verbose",mlmg.cg_verbose );
        pp.query("max_iter",mlmg.max_iter);
        pp.query("max_fmg_iter",mlmg.max_fmg_iter);
        pp.query("agglomeration",mlmg.agglomeration);
        pp.query("consolidation",mlmg.consolidation);
        pp.query("max_coarsening_level",mlmg.max_coarsening_level);
        pp.query("fixed_iter",mlmg.fixed_iter);
    }
    

    // 
    // Initialize geometry and grids
    //    
    Vector<Geometry> geom;
  	Vector<BoxArray> cgrids, ngrids;
 	Vector<DistributionMapping> dmap;
  	Vector<MultiFab> solution, rhs;
 	geom.resize(mesh.nlevels);
 	cgrids.resize(mesh.nlevels);
 	ngrids.resize(mesh.nlevels);
 	dmap.resize(mesh.nlevels);
 	solution.resize(mesh.nlevels);
 	rhs.resize(mesh.nlevels);
	RealBox rb({AMREX_D_DECL(-0.5,-0.5,-0.5)},
	          {AMREX_D_DECL(0.5,0.5,0.5)});
	Geometry::Setup(&rb, 0);
	Box NDomain(IntVect{AMREX_D_DECL(0,0,0)}, 
                IntVect{AMREX_D_DECL(mesh.nnodes,mesh.nnodes,mesh.nnodes)}, 
                IntVect::TheNodeVector());
	Box CDomain = convert(NDomain, IntVect::TheCellVector());

    //
    // Refine the grid
    // 
	Box domain = CDomain;
 	for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
 		{
 			geom[ilev].define(domain);
 			domain.refine(2);
 		}
	Box cdomain = CDomain;
 	for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
	{
		cgrids[ilev].define(cdomain);
		cgrids[ilev].maxSize(mesh.max_grid_size); // TODO
		cdomain.grow(-mesh.nnodes/4); 
		cdomain.refine(2); 
		ngrids[ilev] = cgrids[ilev];
		ngrids[ilev].convert(IntVect::TheNodeVector());
	}

    //
    // Initialize the solution and rhs fabs.
    // Initialize the RHS fab to the function:
    //    RHS[0] = x1*(1-x1) * x2(1-x2) * x3(1-x3)
    //    RHS[1] = 0 
    //    RHS[2] = 0 ... etc
    //
    int nghost = 2;
 	for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
 	{
 		dmap   [ilev].define(cgrids[ilev]);
 		solution[ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost); 
        solution[ilev].setVal(0.0);
        solution[ilev].setMultiGhost(true);
 		rhs     [ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        rhs     [ilev].setVal(0.0);
        rhs     [ilev].setMultiGhost(true);
           
	    Box domain(geom[ilev].Domain());
        const Real AMREX_D_DECL( dx = geom[ilev].CellSize()[0],
                                 dy = geom[ilev].CellSize()[1],
                                 dz = geom[ilev].CellSize()[2]);
        const Real AMREX_D_DECL( minx = geom[ilev].ProbLo()[0],
                                 miny = geom[ilev].ProbLo()[1],
                                 minz = geom[ilev].ProbLo()[2]);
	    domain.convert(IntVect::TheNodeVector());
	    domain.grow(-1); // Shrink domain so we don't operate on any boundaries            
        for (MFIter mfi(solution[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    	{
    		Box bx = mfi.tilebox();
    		bx.grow(1);        // Expand to cover first layer of ghost nodes
    		bx = bx & domain;  // Take intersection of box and the problem domain
		
    		Array4<Real> const& RHS  = rhs[ilev].array(mfi);
    		for (int n = 0; n < op.ncomp; n++)
    			ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    
                    Real AMREX_D_DECL(x1 = i*dx + minx,
                                      x2 = j*dy + miny, 
                                      x3 = k*dz + minz);

                    if (n==0) RHS(i,j,k,n) = AMREX_D_TERM(   (x1-0.5)*(x1+0.5),
                                                           * (x2-0.5)*(x2+0.5),
                                                           * (x3-0.5)*(x3+0.5));
                    else RHS(i,j,k,n) = 0.0;
    			});         
 	    }
        rhs[ilev].FillBoundary(false,true);
    }
         
    // 
    // Set params to be passed to MLMG solver
    //
    LPInfo info;
    if (mlmg.agglomeration >= 0)        info.setAgglomeration(mlmg.agglomeration);
    if (mlmg.consolidation >= 0)        info.setConsolidation(mlmg.consolidation);
    if (mlmg.max_coarsening_level >= 0) info.setMaxCoarseningLevel(mlmg.max_coarsening_level);
    
    //
    // Initialize the MCNodalLinOp linear operator 
    // (see ./MCNodalLinOp.cpp, ./MCNodalLinOp.H for implementation)
    //
    MCNodalLinOp linop;
    linop.setNComp(op.ncomp);
    linop.setCoeff(op.coeff);
    linop.define(geom,cgrids,dmap,info);
    linop.setDomainBC({AMREX_D_DECL(amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet)},
                      {AMREX_D_DECL(amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet)});
    for (int ilev = 0; ilev < mesh.nlevels; ++ilev) linop.setLevelBC(ilev,&solution[ilev]);

    //
    // Initialize the MLMG solver
    //
    MLMG solver(linop);
    if (mlmg.verbose >= 0)     solver.setVerbose(mlmg.verbose);
    if (mlmg.cg_verbose >= 0)  solver.setCGVerbose(mlmg.cg_verbose);
    if (mlmg.fixed_iter >= 0)  solver.setFixedIter(mlmg.fixed_iter);
    if (mlmg.max_iter >= 0)    solver.setMaxIter(mlmg.max_iter);
    if (mlmg.max_fmg_iter >= 0)solver.setMaxFmgIter(mlmg.max_fmg_iter);
    // IMPORTANT! Use the "CFStrategy::ghostnodes" strategy to avoid
    // having to implement a complicated "reflux" routine!
    solver.setCFStrategy(MLMG::CFStrategy::ghostnodes);
    
    //
    // Perform the solve
    //
    Real tol_rel = 1E-8, tol_abs = 1E-8;
    solver.solve(GetVecOfPtrs(solution),GetVecOfConstPtrs(rhs),tol_rel,tol_abs);

    //
    // Write the output to ./solution
    //
    WriteMLMF ("solution",GetVecOfConstPtrs(solution),geom);

    }
    Finalize();
}

