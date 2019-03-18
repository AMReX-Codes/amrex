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

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

using namespace amrex;

int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
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
    amrex::Print() << "**********************************\n"; 
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

    Gpu::setLaunchRegion(true);

    Real cells = 0;
    {
       MultiFab x(ba, dm, Ncomp, Nghost);
       MultiFab y(ba, dm, Ncomp, Nghost);

       Real x_val = 1.0;
       Real y_val = 1.0;

       x.setVal(x_val);
       y.setVal(y_val);

       Real strt_time = amrex::second();
       cells = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
       Real end_time = amrex::second();

       amrex::Print() << "Total number of cells: " << cells << "." << std::endl;
       amrex::Print() << " calculated in " << (end_time-strt_time) << " seconds."<< std::endl << std::endl;
    }

    {
       MultiFab x(ba, dm, Ncomp, Nghost);
       MultiFab y(ba, dm, Ncomp, Nghost);

       Real x_val = 2.0;
       Real y_val = 4.0;

       x.setVal(x_val);
       y.setVal(y_val);

       Real strt_time = amrex::second();
       Real dot_result = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
       Real end_time = amrex::second();
 
       amrex::Print() << "GPU: " << x_val << " dot " << y_val << " = " << dot_result << std::endl;
       amrex::Print() << "Expected value: " << (x_val * y_val * cells) << std::endl;
       amrex::Print() << "Calculated in " << (end_time-strt_time) << " seconds."<< std::endl << std::endl;
    }

    Gpu::setLaunchRegion(false);
    {
       MultiFab x(ba, dm, Ncomp, Nghost);
       MultiFab y(ba, dm, Ncomp, Nghost);

       Real x_val = 2.0;
       Real y_val = 4.0;

       x.setVal(x_val);
       y.setVal(y_val);

       Real strt_time = amrex::second();
       Real dot_result = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
       Real end_time = amrex::second();
 
       amrex::Print() << "CPU: " << x_val << " dot " << y_val << " = " << dot_result << std::endl;
       amrex::Print() << "Expected value: " << (x_val * y_val * cells) << std::endl;
       amrex::Print() << "Calculated in " << (end_time-strt_time) << " seconds."<< std::endl << std::endl;
    }

    Gpu::setLaunchRegion(true);
    {
       MultiFab x(ba, dm, Ncomp, Nghost);
       x.setVal(0.0);

       for (int k = 0; k < 10; k++)
           for (int j = 0; j < 10; j++)
               for (int i = 0; i < 10; i++)
               {
                   IntVect iv(i,j,k);
                   x[0](iv) = 0.01*(i-5) + 0.1*(j-5) + (k-5);
               }

       Real min_val = 0.01*(-5) + 0.1*(-5) + -5;
       Real max_val = 0.01*(4)  + 0.1*(4)  + 4;

       Real strt_time = amrex::second();
       Real min = x.min(0, Nghost); 
       Real max = x.max(0, Nghost); 
       Real end_time = amrex::second();
 
       amrex::Print() << "GPU, expected min/max: " << min_val << "/" << max_val << std::endl;
       amrex::Print() << "Calculatd min/max: " << min << "/" << max << std::endl;
       amrex::Print() << "Calculated in " << (end_time-strt_time) << " seconds."<< std::endl << std::endl;
    }

    Gpu::setLaunchRegion(false);
    {
       MultiFab x(ba, dm, Ncomp, Nghost);
       x.setVal(0.0);

       for (int k = 0; k < 10; k++)
           for (int j = 0; j < 10; j++)
               for (int i = 0; i < 10; i++)
               {
                   IntVect iv(i,j,k);
                   x[0](iv) = 0.01*(i-5) + 0.1*(j-5) + (k-5);
               }

       Real min_val = 0.01*(-5) + 0.1*(-5) + -5;
       Real max_val = 0.01*(4)  + 0.1*(4)  + 4;

       Real strt_time = amrex::second();
       Real min = x.min(0, Nghost); 
       Real max = x.max(0, Nghost); 
       Real end_time = amrex::second();
 
       amrex::Print() << "CPU, expected min/max: " << min_val << "/" << max_val << std::endl;
       amrex::Print() << "Calculatd min/max: " << min << "/" << max << std::endl;
       amrex::Print() << "Calculated in " << (end_time-strt_time) << " seconds."<< std::endl << std::endl;
    }


    amrex::Print() << std::endl << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl << std::endl;
    }

    amrex::Finalize();
}
