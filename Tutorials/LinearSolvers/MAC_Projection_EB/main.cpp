#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MacProjector.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void write_plotfile(const Geometry& geom, const MultiFab& plotmf)
{
    std::stringstream sstream;
    sstream << "plt00000";
    std::string plotfile_name = sstream.str();

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

        amrex::Vector<int> ob_id;

        // read parameters
        {
            ParmParse pp;
            pp.query("mg_verbose", mg_verbose);
            pp.query("cg_verbose", cg_verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("use_hypre", use_hypre);

            pp.queryarr("obstacles", ob_id);

            pp.query("obstacle_radius", obstacle_radius);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1) 
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");

        int n_cell_x = 2*n_cell;
        int n_cell_y =   n_cell;
        int n_cell_z =   n_cell/8;
        int num_obstacles;

        if (ob_id.empty())
        {
           amrex::Print() << " **************************************************** "     << std::endl;
           amrex::Print() << " You didn't specify any obstacles -- please try again " << std::endl;
           amrex::Print() << " ****************************************************\n "     << std::endl;
           exit(0);

        } else {

           num_obstacles = ob_id.size();

           if (num_obstacles > 9)
           {
              amrex::Print() << " **************************************************** "     << std::endl;
              amrex::Print() << " We only have 9 possible obstacles " << std::endl;
              amrex::Print() << " You specified too many -- please try again " << std::endl;
              amrex::Print() << " ****************************************************\n "     << std::endl;
              exit(0);
           } 

           for (int i = 0; i < num_obstacles; i++) 
              if (ob_id[i] < 0 || ob_id[i] > 8)
              {
                 amrex::Print() << " **************************************************** "     << std::endl;
                 amrex::Print() << " The obstacles must be identified using integers from 0 through 8 (inclusive) " << std::endl;
                 amrex::Print() << " You specified an invalid obstacle -- please try again " << std::endl;
                 amrex::Print() << " ****************************************************\n "     << std::endl;
                 exit(0);
              }

           amrex::Print() << " \n********************************************************************" << std::endl; 
           amrex::Print() << " You specified " << num_obstacles << " objects in the domain: ";
              for (int i = 0; i < num_obstacles; i++) 
                  amrex::Print() << ob_id[i] << " ";
             amrex::Print() << std::endl;
           amrex::Print() << " ********************************************************************" << std::endl; 
        } 

        Real zlen = 0.125;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(2.0,1.0,zlen)});

            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;
        MultiFab plotfile_mf;

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
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 0], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 1], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 2], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 3], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 4], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 5], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 6], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 7], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 8], false)};

        switch(num_obstacles) {

           case 1:
              {
              auto gshop1 = EB2::makeShop(obstacles[ob_id[0]]);
              EB2::Build(gshop1, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 2:
              {
              amrex::Print() << "Objects " << ob_id[0] << " " << ob_id[1] << std::endl;
              auto all2 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]]);
              auto gshop2  = EB2::makeShop(all2);
              EB2::Build(gshop2, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 3:
              {
              auto all3 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto gshop3  = EB2::makeShop(all3);
              EB2::Build(gshop3, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 4:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto all     = EB2::makeUnion(group_1,obstacles[ob_id[3]]);
              auto gshop4  = EB2::makeShop(all);
              EB2::Build(gshop4, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 5:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]]);
              auto all     = EB2::makeUnion(group_1,group_2);
              auto gshop5  = EB2::makeShop(all);
              EB2::Build(gshop5, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 6:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto all     = EB2::makeUnion(group_1,group_2);
              auto gshop6  = EB2::makeShop(all);
              EB2::Build(gshop6, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 7:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = obstacles[ob_id[6]];
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop7  = EB2::makeShop(all);
              EB2::Build(gshop7, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 8:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = EB2::makeUnion(obstacles[ob_id[6]],obstacles[ob_id[7]]);
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop8  = EB2::makeShop(all);
              EB2::Build(gshop8, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 9:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = EB2::makeUnion(obstacles[ob_id[6]],obstacles[ob_id[7]],obstacles[ob_id[8]]);
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop9  = EB2::makeShop(all);
              EB2::Build(gshop9, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           default:;
        }
   
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);
   
        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;
  
        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};
 
        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

	// allocate face-centered velocities and face-centered beta coefficient
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1,
			      MFInfo(), factory);
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0,
			      MFInfo(), factory);
            beta[idim].setVal(1.0);  // set beta to 1
        }

	// If we want to supply a non-zero S we must allocate and fill it outside the solver
        // MultiFab S(grids, dmap, 1, 0, MFInfo(), factory);
	// Set S here ... 

        // store plotfile variables; velocity and processor id
        plotfile_mf.define(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), factory);

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
        if (use_hypre) 
            lp_info.setMaxCoarseningLevel(0);

        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                             {amrex::GetArrOfConstPtrs(beta)}, // beta
                             {geom},                           // the geometry object
                             lp_info);                         // structure for passing info to the operator
	  
	// Here we specifiy the desired divergence S
	// MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // face-based velocity
	//                      {amrex::GetArrOfConstPtrs(beta)}, // beta
	//                      {geom},                           // the geometry object
	//                      lp_info,                          // structure for passing info to the operator
	//                      {&S});                            // defines the specified RHS divergence

        // Set bottom-solver to use hypre instead of native BiCGStab 
        if (use_hypre) 
           macproj.setBottomSolver(MLMG::BottomSolver::hypre);
	
	// Hard-wire the boundary conditions to be Neumann on the low x-face, Dirichlet
	// on the high x-face, and periodic in the other two directions  
	// (the first argument is for the low end, the second is for the high end)
        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Periodic,
                                          LinOpBCType::Periodic)},
	                    {AMREX_D_DECL(LinOpBCType::Dirichlet,
					  LinOpBCType::Periodic,
					  LinOpBCType::Periodic)});

        macproj.setVerbose(mg_verbose);
        macproj.setCGVerbose(cg_verbose);
	
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

	// Solve for phi and subtract from the velocity to make it divergence-free
        macproj.project(reltol,abstol);
	
	// If we want to use phi elsewhere, we can pass in an array in which to return the solution
	// MultiFab phi_inout(grids, dmap, 1, 1, MFInfo(), factory);	
	// macproj.project({&phi_inout},reltol,abstol);

        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " Done!" << std::endl;
        amrex::Print() << "******************************************************************** \n" << std::endl; 

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }

        // copy processor id into plotfile_mf
	plotfile_mf.setVal(ParallelDescriptor::MyProc(), 0, 1);
	
        // copy velocity into plotfile
        average_face_to_cellcenter(plotfile_mf,1,amrex::GetArrOfConstPtrs(vel));

        write_plotfile(geom, plotfile_mf); 
    }
  
    Real stop_time = amrex::second() - strt_time;
    amrex::Print() << "Total run time " << stop_time << std::endl;

    amrex::Finalize();
}
