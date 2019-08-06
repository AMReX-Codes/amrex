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
    if (argc < 2) 
    {
        std::cout << "Missing input file" << std::endl;
        exit(-1);
    }
    Initialize(argc, argv);
    
    struct {
        int nlevels = 3;
        int nnodes = 32;
    } mesh;
    {
        ParmParse pp("mesh");
        pp.query("nlevels",mesh.nlevels);
        pp.query("nnodes",mesh.nnodes);
    }
    
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
    
    struct {
        int ncomp=1;
        Vector<Real> coeff = {1.0};
    } op;
    {
        ParmParse pp("op");
        pp.query("ncomp",op.ncomp);
        pp.queryarr("coeff",op.coeff);
    }

    
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

	RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)},
			          {AMREX_D_DECL(1.0,1.0,1.0)});
	Geometry::Setup(&rb, 0);

	Box NDomain(IntVect{AMREX_D_DECL(0,0,0)}, 
                IntVect{AMREX_D_DECL(mesh.nnodes,mesh.nnodes,mesh.nnodes)}, 
                IntVect::TheNodeVector());
	Box CDomain = convert(NDomain, IntVect::TheCellVector());

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
		cgrids[ilev].maxSize(10000); // TODO
		cdomain.grow(-mesh.nnodes/4); 
		cdomain.refine(2); 
		ngrids[ilev] = cgrids[ilev];
		ngrids[ilev].convert(IntVect::TheNodeVector());
	}

    int nghost = 2;
 	for (int ilev = 0; ilev < mesh.nlevels; ++ilev)
 	{
 		dmap   [ilev].define(cgrids[ilev]);
 		solution[ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost); 
        solution[ilev].setVal(0.0);
 		rhs     [ilev].define(ngrids[ilev], dmap[ilev], op.ncomp, nghost);
        rhs     [ilev].setVal(0.0);
           
	    Box domain(geom[ilev].Domain());
        const Real* DX = geom[ilev].CellSize();
	    domain.convert(IntVect::TheNodeVector());
	    domain.grow(-1); // Shrink domain so we don't operate on any boundaries            
        for (MFIter mfi(solution[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    	{
    		Box bx = mfi.tilebox();
    		bx.grow(1);        // Expand to cover first layer of ghost nodes
    		bx = bx & domain;  // Take intersection of box and the problem domain
		
	   		Array4<Real> const& SOL  = solution[ilev].array(mfi);
    		Array4<Real> const& RHS  = rhs[ilev].array(mfi);
    		for (int n = 0; n < op.ncomp; n++)
    			ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    
                    Real x1 = i*DX[0], x2 = j*DX[1], x3 = k*DX[2];

                    if (n==0) RHS(i,j,k,n) = x1*(1.0 - x1) * x2 * (1.0 - x2) * x3 * (1.0 - x3);
                    else RHS(i,j,k,n) = 0.0;
    			});         
 	    }
    }
         

    LPInfo info;
    if (mlmg.agglomeration >= 0)        info.setAgglomeration(mlmg.agglomeration);
    if (mlmg.consolidation >= 0)        info.setConsolidation(mlmg.consolidation);
    if (mlmg.max_coarsening_level >= 0) info.setMaxCoarseningLevel(mlmg.max_coarsening_level);
    
    MCNodalLinOp linop;
    linop.setNComp(op.ncomp);
    linop.setCoeff(op.coeff);
    linop.define(geom,cgrids,dmap,info);
//    linop.setDomainBC({amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet},
//                      {amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet,amrex::MLLinOp::BCType::Dirichlet});
//    for (int ilev = 0; ilev < mesh.nlevels; ++ilev) linop.setLevelBC(ilev,&solution[ilev]);

    MLMG solver(linop);
    solver.setCFStrategy(MLMG::CFStrategy::ghostnodes);
    //solver.setBottomSolver(amrex::BottomSolver::smoother);
    if (mlmg.verbose >= 0)     solver.setVerbose(mlmg.verbose);
    if (mlmg.cg_verbose >= 0)  solver.setCGVerbose(mlmg.cg_verbose);
    if (mlmg.fixed_iter >= 0)  solver.setFixedIter(mlmg.fixed_iter);
    if (mlmg.max_iter >= 0)    solver.setMaxIter(mlmg.max_iter);
    if (mlmg.max_fmg_iter >= 0)solver.setMaxFmgIter(mlmg.max_fmg_iter);
    

    Real tol_rel = 1E-8, tol_abs = 1E-8;
    solver.solve(GetVecOfPtrs(solution),GetVecOfConstPtrs(rhs),tol_rel,tol_abs);

    WriteMLMF ("solution",GetVecOfConstPtrs(solution),geom);
    WriteMLMF ("rhs",GetVecOfConstPtrs(rhs),geom);
    
}

