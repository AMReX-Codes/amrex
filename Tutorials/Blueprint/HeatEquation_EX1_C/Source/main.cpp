//
// Demonstrates how to generate of a Conduit Mesh Blueprint 
// description of an AMReX Single Level dataset and render this
// data in situ using with ALPINE Ascent.
// 

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Conduit_Blueprint.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "myfunc.H"
#include "myfunc_F.H"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>

#include <ascent.hpp>

using namespace conduit;
using namespace ascent;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of 
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        pp.queryarr("is_periodic", is_periodic);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 1;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi());
    }

    // compute the time step
    const Real* dx = geom.CellSize();
    Real dt = 0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);

    // time = starting time in the simulation
    Real time = 0.0;

    int rank   = ParallelDescriptor::MyProc();  // Return the rank
    int ntasks = ParallelDescriptor::NProcs();  // Return the number of processes
    
    /////////////////////////////
    // Setup Ascent
    /////////////////////////////
    // Create an instance of Ascent
    Ascent ascent;
    Node open_opts;

    // for the MPI case, provide the mpi comm
#ifdef BL_USE_MPI
    open_opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

    ascent.open(open_opts);

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)

    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);

        ///////////////////////////////////////////////////////////////////
        // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
        ///////////////////////////////////////////////////////////////////
        conduit::Node bp_mesh;
        SingleLevelToBlueprint( phi_new, {"phi"}, geom, time, n, bp_mesh);

        
        ///////////////////////////////////////////////////////////////////
        // Save the Blueprint Mesh to a set of files that we can 
        // view in VisIt. 
        // (For debugging and to demonstrate how to do this w/o Ascent)
        ///////////////////////////////////////////////////////////////////
        WriteBlueprintFiles( bp_mesh, "bp_heateq_ex_", n);

        ///////////////////////////////////////////////////////////////////
        // Render with Ascent
        ///////////////////////////////////////////////////////////////////

        Node scenes;
        scenes["s1/plots/p1/type"] = "pseudocolor";
        scenes["s1/plots/p1/field"] = "phi";
        // Set the output file name (ascent will add ".png")
        const std::string& png_out = amrex::Concatenate("ascent_render_",0,5);
        scenes["s1/image_prefix"] = png_out;

        // setup actions
        Node actions;
        Node &add_act = actions.append();
        add_act["action"] = "add_scenes";
        add_act["scenes"] = scenes;
        actions.append()["action"] = "execute";
        actions.append()["action"] = "reset";


        ascent.publish(bp_mesh);
        ascent.execute(actions);
    }

    // build the flux multifabs
    std::array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 0);
    }

    

    
    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);

        // new_phi = old_phi + dt * (something)
        advance(phi_old, phi_new, flux, dt, geom); 
        time = time + dt;
        
        
        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data and write a conduit blueprint
        // file as well
        if ( plot_int >= 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);

            ///////////////////////////////////////////////////////////////////
            // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
            ///////////////////////////////////////////////////////////////////
            conduit::Node bp_mesh;
            SingleLevelToBlueprint( phi_new, {"phi"}, geom, time, n, bp_mesh);

            
            ///////////////////////////////////////////////////////////////////
            // Save the Blueprint Mesh to a set of files that we can 
            // view in VisIt. 
            // (For debugging and to demonstrate how to do this w/o Ascent)
            ///////////////////////////////////////////////////////////////////
            WriteBlueprintFiles( bp_mesh, "bp_heateq_ex_", n);

            ///////////////////////////////////////////////////////////////////
            // Render with Ascent
            ///////////////////////////////////////////////////////////////////
            
            // add a scene with a pseudocolor plot
            Node scenes;
            scenes["s1/plots/p1/type"] = "pseudocolor";
            scenes["s1/plots/p1/field"] = "phi";
            // Set the output file name (ascent will add ".png")
            const std::string& png_out = amrex::Concatenate("ascent_render_",n,5);
            scenes["s1/image_prefix"] = png_out;

            // setup actions
            Node actions;
            Node &add_act = actions.append();
            add_act["action"] = "add_scenes";
            add_act["scenes"] = scenes;
            actions.append()["action"] = "execute";
            actions.append()["action"] = "reset";


            ascent.publish(bp_mesh);
            ascent.execute(actions);

        }
    }

    ascent.close();

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
