#include <cuda_runtime.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
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
    int n_cell, max_grid_size;
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

    {
        BL_PROFILE("Standard");
        MultiFab x(ba, dm, Ncomp, Nghost);

        for (MFIter mfi(x); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.validbox();
            Array4<Real> a = x.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k) = 1.0;
            });
        }

        amrex::Print() << "Sum = " << x.sum() << std::endl;
    }

    {
        BL_PROFILE("cudaGraph-iter");
        MultiFab x(ba, dm, Ncomp, Nghost);

        cudaGraph_t     graph[x.local_size()];
        cudaGraphExec_t graphExec[x.local_size()];

        for (MFIter mfi(x); mfi.isValid(); ++mfi)
        {
            AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Cuda::Device::cudaStream()));

            const Box bx = mfi.validbox();
            Array4<Real> a = x.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k) = 2.0;
            });

            AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
        }

        for (int i = 0; i<x.local_size(); ++i)
        {
            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
        }

        for (MFIter mfi(x); mfi.isValid(); ++mfi)
        {
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[mfi.LocalIndex()], amrex::Cuda::Device::cudaStream())); 
        }

        amrex::Gpu::Device::synchronize();

        amrex::Print() << "Sum = " << x.sum() << std::endl;
    }

    {
        BL_PROFILE("cudaGraph-stream");
        MultiFab x(ba, dm, Ncomp, Nghost);

        cudaGraph_t     graph[amrex::Gpu::Device::numCudaStreams()];
        cudaGraphExec_t graphExec[amrex::Gpu::Device::numCudaStreams()];

        for (MFIter mfi(x); mfi.isValid(); ++mfi)
        {
            if (mfi.LocalIndex() == 0)
            {
                for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
                {
                    amrex::Gpu::Device::setStreamIndex(i);
                    AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Cuda::Device::cudaStream()));
                }
                amrex::Gpu::Device::setStreamIndex(mfi.tileIndex());
            } 

            const Box bx = mfi.validbox();
            Array4<Real> a = x.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k) = 3.0;
            });

            if (mfi.LocalIndex() == (x.local_size() - 1) )
            { 
                for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
                {
                    amrex::Gpu::Device::setStreamIndex(i); 
                    AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[i])));
                }
            }
        }

        for (int i = 0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
        {
            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
        }

        for (int i = 0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
        {
            amrex::Gpu::Device::setStreamIndex(i); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[i], amrex::Cuda::Device::cudaStream())); 
        }

        amrex::Gpu::Device::resetStreamIndex();

        amrex::Gpu::Device::synchronize();

        amrex::Print() << "Sum = " << x.sum() << std::endl;
    }

    amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
