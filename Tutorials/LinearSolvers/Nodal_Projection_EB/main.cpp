#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLMG.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void write_plotfile(const Geometry& geom, const MultiFab& plotmf)
{
    std::string plotfile_name("plt00000");

    amrex::Print() << "Writing " << plotfile_name << std::endl;    
    
#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc" ,"xvel", "yvel" },
                                     geom, 0.0, 0);
#elif (AMREX_SPACEDIM == 3)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc", "xvel", "yvel", "zvel" },
                                     geom, 0.0, 0);
#endif
}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    Real strt_time = amrex::second();
    
    {
        int mg_verbose = 0;
        int cg_verbose = 0;
        int n_cell = 128;
        int max_grid_size = 32;
        int use_hypre  = 0;

        Real obstacle_radius = 0.10;

        // read parameters
        {
            ParmParse pp;
            pp.query("mg_verbose", mg_verbose);
            pp.query("cg_verbose", cg_verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("use_hypre", use_hypre);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1) 
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");

        int n_cell_y =   n_cell;
        int n_cell_x = 2*n_cell;
        int n_cell_z =   n_cell/8;

        Real ylen = 1.0;
        Real xlen = 2.0 * ylen;
        Real zlen = ylen / 8.0;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(xlen,ylen,zlen)});
            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain, rb, CoordSys::cartesian, isp);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        amrex::Vector<amrex::RealArray> obstacle_center = {
            {AMREX_D_DECL(0.3,0.2,0.5)},
            {AMREX_D_DECL(0.3,0.5,0.5)},
            {AMREX_D_DECL(0.3,0.8,0.5)},
            {AMREX_D_DECL(0.7,0.25,0.5)},
            {AMREX_D_DECL(0.7,0.60,0.5)},
            {AMREX_D_DECL(0.7,0.85,0.5)},
            {AMREX_D_DECL(1.1,0.2,0.5)},
            {AMREX_D_DECL(1.1,0.5,0.5)},
            {AMREX_D_DECL(1.1,0.8,0.5)}};

        int direction =  2;
        Real height   = -1.0;  // Putting a negative number for height means it extends beyond the domain

        // The "false" below is the boolean that determines if the fluid is inside ("true") or 
        //     outside ("false") the object(s)

        Array<EB2::CylinderIF,9> obstacles{
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 0], false),
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 1], false),
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 2], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 3], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 4], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 5], false),
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 6], false),
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 7], false),
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 8], false)};

        auto group_1 = EB2::makeUnion(obstacles[0],obstacles[1],obstacles[2]);
        auto group_2 = EB2::makeUnion(obstacles[3],obstacles[4],obstacles[5]);
        auto group_3 = EB2::makeUnion(obstacles[6],obstacles[7],obstacles[8]);
        auto all     = EB2::makeUnion(group_1,group_2,group_3);
        auto gshop9  = EB2::makeShop(all);
        EB2::Build(gshop9, geom, required_coarsening_level, max_coarsening_level);
   
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);
   
        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;
  
        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};
 
        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

        //
        // Given a cell-centered velocity (vel) field, a cell-centered
        // scalar field (sigma) field, and a source term S (either node-
        // or cell-centered) solve:
        //
        //   div( sigma * grad(phi) ) = div(vel) - S
        //
        // and then perform the projection:
        //
        //     vel = vel - sigma * grad(phi)
        // 

        //
        //  Create the cell-centered velocity field we want to project  
        //
        MultiFab vel(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);

        // Set velocity field to (1,0,0) including ghost cells for this example
        vel.setVal(1.0, 0, 1, 1);
        vel.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);

        //
        // Setup linear operator, AKA the nodal Laplacian
        // 
        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
        // if (use_hypre_as_full_solver)
        //    lp_info.setMaxCoarseningLevel(0);

        MLNodeLaplacian matrix({geom}, {grids}, {dmap}, lp_info,
                               Vector<EBFArrayBoxFactory const*>{&factory});

        // Set boundary conditions.
        // Here we use Neumann on the low x-face, Dirichlet on the high x-face,
        // and periodic in the other two directions
        // (the first argument is for the low end, the second is for the high end)
        // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)
        matrix.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)},
                           {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                          LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)});

        // Set matrix attributes to be used by MLMG solver
        matrix.setGaussSeidel(true);
        matrix.setHarmonicAverage(false);

        //
        // Compute RHS 
        //
        // NOTE: it's up to the user to compute the RHS. as opposed to the MAC projection 
        //
        // NOTE: do this operation AFTER setting up the linear operator so
        //       that compRHS method can be used
        // 

        // RHS is nodal
        const BoxArray & nd_grids = amrex::convert(grids, IntVect::TheNodeVector()); // nodal grids
 
        // Multifab to host RHS
        MultiFab rhs(nd_grids, dmap, 1, 1, MFInfo(), factory);
 
        // Cell-centered contributions to RHS
        MultiFab S_cc(grids, dmap, 1, 1, MFInfo(), factory);
        S_cc.setVal(0.0); // Set it to zero for this example
  
       // Node-centered contributions to RHS
        MultiFab S_nd(nd_grids, dmap, 1, 1, MFInfo(), factory);
        S_nd.setVal(0.0); // Set it to zero for this example 
  
        // Compute RHS -- vel must be cell-centered
        matrix.compRHS({&rhs}, {&vel}, {&S_nd}, {&S_cc});
 
        //
        // Create the cell-centered sigma field and set it to 1 for this example
        //
        MultiFab sigma(grids, dmap, 1, 1, MFInfo(), factory);
        sigma.setVal(1.0);
  
        // Set sigma 
        matrix.setSigma(0, sigma);
  
        //
        // Create node-centered phi
        //
        MultiFab phi(nd_grids, dmap, 1, 1, MFInfo(), factory);
        phi.setVal(0.0);
 
        //
        // Setup MLMG solver
        //
        MLMG nodal_solver(matrix);
 
        // We can specify the maximum number of iterations
        // nodal_solver.setMaxIter(nodal_mg_maxiter);
        // nodal_solver.setCGMaxIter(nodal_mg_cg_maxiter);
 
        nodal_solver.setVerbose(mg_verbose);
        nodal_solver.setCGVerbose(cg_verbose);

        // Set bottom-solver to use hypre instead of native BiCGStab 
        //   ( we could also have set this to cg, bicgcg, cgbicg)
        // if (use_hypre_as_full_solver || use_hypre_as_bottom_solver) 
        //    nodal_solver.setBottomSolver(MLMG::BottomSolver::hypre);
 
        // Define the relative tolerance
        Real reltol = 1.e-8;
 
        // Define the absolute tolerance; note that this argument is optional
        Real abstol = 1.e-15;
 	
        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " Let's project the initial velocity to find " << std::endl;
        amrex::Print() << "   the flow field around the obstacles ... " << std::endl;
        amrex::Print() << " The domain has " << n_cell_x << " cells in the x-direction "          << std::endl;
        amrex::Print() << " The maximum grid size is " << max_grid_size                             << std::endl;  
        amrex::Print() << "******************************************************************** \n" << std::endl; 
 
        //
        // Solve div( sigma * grad(phi) ) = RHS
        //
        nodal_solver.solve( {&phi}, {&rhs}, reltol, abstol);
 
        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " Done solving the equation " << std::endl;
        amrex::Print() << " ... now subtracting off sigmna grad phi from vel" << std::endl;
 
        //
        // Create cell-centered multifab to hold value of -sigma*grad(phi) at cell-centers
        // 
        MultiFab fluxes(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
        fluxes.setVal(0.0);
 
        // Get fluxes from solver
        nodal_solver.getFluxes( {&fluxes} );
 
        //
        // Apply projection explicitly --  vel = vel - sigma * grad(phi)  
        // 
        MultiFab::Add( vel, fluxes, 0, 0, AMREX_SPACEDIM, 0);
 
        amrex::Print() << " ... now done with full projection operation" << std::endl;
        amrex::Print() << "******************************************************************** \n" << std::endl; 
 
        // Store plotfile variables; velocity and processor id
        MultiFab plotfile_mf(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), factory);
 
        // copy processor id into plotfile_mf
        plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);
        plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);
 
        // copy velocity into plotfile
        MultiFab::Copy(plotfile_mf, vel, 0, 1, AMREX_SPACEDIM, 0);

        write_plotfile(geom, plotfile_mf); 
    }
  
    Real stop_time = amrex::second() - strt_time;
    amrex::Print() << "Total run time " << stop_time << std::endl;

    amrex::Finalize();
}
