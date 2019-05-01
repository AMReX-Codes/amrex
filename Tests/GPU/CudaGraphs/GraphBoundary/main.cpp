#include <cuda_device_runtime_api.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

#include <Prob.H>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
    amrex::Print() << "amrex::Initialize complete." << "\n";

    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices = 0;
#ifdef AMREX_USE_CUDA
    cudaGetDeviceCount(&devices);
#endif
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    amrex::Print() << "**********************************\n\n"; 
    // ===================================

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

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
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // Periodic in all directions by default
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
 
    amrex::Print() << std::endl;
    amrex::Print() << " Domain size: " << n_cell << " cells." << std::endl;
    amrex::Print() << " Max grid size: " << max_grid_size << " cells." << std::endl;
    amrex::Print() << " Boxes: " << ba.size() << std::endl << std::endl;
    amrex::Print() << "========================================================" << std::endl << std::endl;

    MultiFab mf_graph(ba, dm, Ncomp, Nghost);
    MultiFab mf_gpu  (ba, dm, Ncomp, Nghost);
    MultiFab mf_cpu  (ba, dm, Ncomp, Nghost);

    Real start_time, end_time;
//    Real cpu_avg, gpu_avg, graph_avg, graph_init;
    Gpu::setLaunchRegion(true); 
    Gpu::setGraphRegion(true); 


    // With GPUs and Graphs
    {
        // Setup the data outside the FillBoundary timers.
        // Ensures data is moved to GPU.
        for (MFIter mfi(mf_graph); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.validbox();
            Array4<Real> phi = mf_graph[mfi].array(); 
            GeometryData geomData = geom.data();

            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                initdata(tbx, phi, geomData);
            });
         }    

        // First FillBoundary will create graph and run.
        // Timed separately.
        BL_PROFILE_VAR("GRAPH: Create Graph and Run", makeandrungraph);

        start_time = amrex::second();
        mf_graph.FillBoundary(geom.periodicity());
        end_time = amrex::second();

        BL_PROFILE_VAR_STOP(makeandrungraph);

        std::cout << ParallelDescriptor::MyProc() << " : Time for initial FillBoundary (Create & Run Graph) = " << end_time - start_time << std::endl;

        // Run the remainder of the FillBoundarys (nsteps-1)
        BL_PROFILE_VAR("GRAPH: Run Graph", rungraph);

        start_time = amrex::second();
        for (int i=1; i<nsteps; ++i)
        {
            mf_graph.FillBoundary(geom.periodicity());
        }
        end_time = amrex::second();

        BL_PROFILE_VAR_STOP(rungraph);

        std::cout << ParallelDescriptor::MyProc() << " : Time for " << nsteps-1 << " subsequent FillBoundary(s) = " << end_time - start_time << std::endl;
        std::cout << ParallelDescriptor::MyProc() << " : Average time per graph-only FillBoundary = " << (end_time - start_time)/(nsteps) << std::endl;
        std::cout << "========================================================" << std::endl << std::endl;
    }

    Gpu::setLaunchRegion(true); 
    Gpu::setGraphRegion(false); 
    // With GPU and no graphs
    {
        // Setup the data outside the FillBoundary timers.
        // Ensures data is moved to GPU.
        for (MFIter mfi(mf_gpu); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.validbox();
            Array4<Real> phi = mf_gpu[mfi].array(); 
            GeometryData geomData = geom.data();

            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                initdata(tbx, phi, geomData);
            });
         }    

        // Run the nstep FillBoundaries
        BL_PROFILE_VAR("GPU: FillBoundary", GPUFB);

        start_time = amrex::second();
        for (int i=0; i<nsteps; ++i)
        {
            mf_gpu.FillBoundary(geom.periodicity());
        }
        end_time = amrex::second();

        BL_PROFILE_VAR_STOP(GPUFB);

        std::cout << ParallelDescriptor::MyProc() << " : Time for " << nsteps << " FillBoundary(s) = " << end_time - start_time << std::endl;
        std::cout << ParallelDescriptor::MyProc() << " : Average time per GPU FillBoundary = " << (end_time - start_time)/(nsteps) << std::endl;
        std::cout << "========================================================" << std::endl << std::endl;
    }

    Gpu::setLaunchRegion(false); 
    // With CPU
    {
        // Setup the data outside the FillBoundary timers.
        // Ensures data is moved to GPU.
        for (MFIter mfi(mf_cpu); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.validbox();
            Array4<Real> phi = mf_cpu[mfi].array(); 
            GeometryData geomData = geom.data();

            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                initdata(tbx, phi, geomData);
            });
         }    

        // Run the remainder of the FillBoundary's (nsteps-1)
        BL_PROFILE_VAR("CPU: FillBoundary", CPUFB);

        start_time = amrex::second();
        for (int i=0; i<nsteps; ++i)
        {
            mf_cpu.FillBoundary(geom.periodicity());
        }
        end_time = amrex::second();

        BL_PROFILE_VAR_STOP(CPUFB);

        std::cout << ParallelDescriptor::MyProc() << " : Time for " << nsteps << " FillBoundary(s) = " << end_time - start_time << std::endl;
        std::cout << ParallelDescriptor::MyProc() << " : Average time per CPU FillBoundary = " << (end_time - start_time)/(nsteps) << std::endl;
        std::cout << "========================================================" << std::endl << std::endl;
    }

    {
        amrex::Real max_error;
        MultiFab mf_error (ba, dm, Ncomp, Nghost);

        MultiFab::Copy(mf_error, mf_cpu, 0, 0, 1, 1);
        MultiFab::Subtract(mf_error, mf_gpu, 0, 0, 1, 1);
        max_error = mf_error.norm0(0, 1);
        amrex::Print() << "Max difference between CPU and GPU: " << max_error << std::endl; 

        MultiFab::Copy(mf_error, mf_graph, 0, 0, 1, 1);
        MultiFab::Subtract(mf_error, mf_cpu, 0, 0, 1, 1);
        max_error = mf_error.norm0(0, 1);
        amrex::Print() << "Max difference between CPU and Graph: " << max_error << std::endl; 

        std::cout << "========================================================" << std::endl << std::endl;
    }

    // Clean up output and timers.
    // Add "Average FillBoundaries to break even" math.
    // SegFault? ( 256^3 & max_grid_size = 8 )
    // CUDA 10 proofing.  Device Calls, FillBoundary, Graph Definitions.
    // Flag for graph per stream & graph per iter.
    // Parallel Copy?
    // Clean up Graph methodology.

    amrex::Print() << std::endl << "************************************************" << std::endl;
    }

    amrex::Finalize();
}
