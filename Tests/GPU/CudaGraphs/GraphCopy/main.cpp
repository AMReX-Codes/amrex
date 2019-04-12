#include <cuda_runtime.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

using namespace amrex;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

__global__
void copy (amrex::Dim3 lo, amrex::Dim3 len, int ncells,
           amrex::Dim3 offset, amrex::Array4<Real> src, amrex::Array4<Real> dst,
           int scomp, int dcomp, int ncomp)
{

    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        for (int n = 0; n < ncomp; ++n) {
            dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
        }
    }

}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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

        // Malloc value for setval testing.
        Real* val;
        cudaMallocManaged(&val, sizeof(Real));

        // Create the MultiFab and touch the data.
        // Ensures the data in on the GPU for all further testing.
//        MultiFab x(ba, dm, Ncomp, Nghost);
//        MultiFab y(ba, dm, Ncomp, Nghost);
//        x.setVal(1.0);
//        y.setVal(2.0);

        Real points = ba.numPts();

        amrex::Print() << "Testing on " << n_cell << "^3 boxes with max grid size " << max_grid_size << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs

        {
            // Create the MultiFab and touch the data.
            // Ensures the data in on the GPU for all further testing.
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            x.setVal(4.5);
            y.setVal(5.3);

            BL_PROFILE("First Copy");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real> src = x.array(mfi);
                Array4<Real> dst = y.array(mfi);
                Dim3 offset {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
                {
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });
            }

            amrex::Print() << "First sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

// ---------------------------------------

        {
            // Create the MultiFab and touch the data.
            // Ensures the data in on the GPU for all further testing.
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            x.setVal(4.0);
            y.setVal(5.0);

            BL_PROFILE("Lambda Copy");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real> src = x.array(mfi);
                Array4<Real> dst = y.array(mfi);
                Dim3 offset {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
                {
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });
            }

            amrex::Print() << "No Graph Lambda sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

// ---------------------------------------

        {
            // Create the MultiFab and touch the data.
            // Ensures the data in on the GPU for all further testing.
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            x.setVal(1.0);
            y.setVal(2.0);

            BL_PROFILE("Function Copy");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const auto ec = Cuda::ExecutionConfig(ncells);
                const Dim3 offset = {0,0,0};

                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, a, b, 0, 0, 1);
            }

            amrex::Print() << "No Graph Function sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }




// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        cudaGraphExec_t graphExec;

        {
            // Create the MultiFab and touch the data.
            // Ensures the data in on the GPU for all further testing.
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            x.setVal(1e-5);
            y.setVal(0.0);

            BL_PROFILE("cudaGraph");

// --------- Capture each stream in the MFIter loop ----------

            BL_PROFILE_VAR("cudaGraph-create", cgc);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    amrex::Gpu::Device::startGraphStreamRecording();
                } 

                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const auto ec = Cuda::ExecutionConfig(ncells);
                const Dim3 offset = {0,0,0};

                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, a, b, 0, 0, 1); 

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                {
                    graphExec = amrex::Gpu::Device::stopGraphStreamRecording(); 
                }
            }

            BL_PROFILE_VAR_STOP(cgc);

// --------- Launch the graph  ----------

            BL_PROFILE_VAR("cudaGraph-launch", cgl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Graphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph with a different result  ----------

            x.setVal(8.67530e-9);
            y.setVal(0.0);

            BL_PROFILE_VAR("cudaGraph-relaunch", cgrl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgrl);

            amrex::Print() << "Regraphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

/*      Written out version, for reference.
        Above version uses Device::() functions to allow pattern to be placed elsewhere.
        {
            BL_PROFILE("cudaGraph");

            // Create the MultiFab and touch the data.
            // Ensures the data in on the GPU for all further testing.
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            x.setVal(1e-5);
            y.setVal(0.0);

            cudaGraph_t     graph[amrex::Gpu::Device::numCudaStreams()];
//            cudaGraphExec_t graphExec;

// --------- Capture each stream in the MFIter loop ----------

            BL_PROFILE_VAR("cudaGraph-create", cgc);

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
                Array4<Real> b = y.array(mfi);

                int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const auto ec = Cuda::ExecutionConfig(ncells);
                const Dim3 offset = {0,0,0};

                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, a, b, 0, 0, 1); 

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[i])));
                    }
                }
            }

// --------- Build full graph ----------
// --------- Note: Current version assumes there is enough work to put something
// --------- in every stream. Likely needs modification to test and put together
// --------- a final graph only with graphs containing captured data. 

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc);

// --------- Instantiate the graph  ----------

            BL_PROFILE_VAR("cudaGraph-instantiate", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);

// --------- Launch the graph  ----------

            BL_PROFILE_VAR("cudaGraph-launch", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Gpu::Device::resetStreamIndex();

            amrex::Print() << "Graphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph with a different result  ----------

            BL_PROFILE_VAR("cudaGraph-relaunch", cgrl);

            x.setVal(8.67530e-9);
            y.setVal(0.0);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();

            BL_PROFILE_VAR_STOP(cgrl);

            amrex::Gpu::Device::resetStreamIndex();

            amrex::Print() << "Regraphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }
*/
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Screw with memory a bit and show that the graphExec doesn't do what we want.
//      Change this example later to test stored pointers to graph input data.

        {
            MultiFab v(ba, dm, Ncomp, Nghost);
            MultiFab w(ba, dm, Ncomp, Nghost);
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);
            v.setVal(0.0);
            w.setVal(0.0);
            x.setVal(0.15);
            y.setVal(0.25);

            BL_PROFILE("cudaGraph-newMFabs");

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            amrex::Print() << "New MultiFabs = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

        }

        amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
