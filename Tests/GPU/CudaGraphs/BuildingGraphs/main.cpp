#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// MFIterLoop
// Written as a seperate function for easy changes/testing.
void MFIterLoopFunc(const Box &bx, double* val, Array4<Real> &a)
{
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        a(i,j,k) = *val;
    });
}


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
    amrex::Initialize(argc, argv);

    {
        amrex::Print() << "amrex::Initialize complete." << "\n";

        // AMREX_SPACEDIM: number of dimensions
        int n_cell, max_grid_size, Nghost, Ncomp;
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

            // The domain is broken into boxes of size max_grid_size
            Nghost = 0;
            pp.query("nghost", Nghost);

            Ncomp = 1;
            pp.query("ncomp", Ncomp);
        }

        // make BoxArray and Geometry
        BoxArray ba;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "ba" from the single box "bx"
            ba.define(domain);
            // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
            ba.maxSize(max_grid_size);
        }
  
        // How Boxes are distrubuted among MPI processes
        DistributionMapping dm(ba);

        // Malloc value for setval testing.
        Real* val;
        cudaMallocManaged(&val, sizeof(Real));
        *val = 0.0;

        // Create the MultiFab and touch the data.
        // Ensures the data in on the GPU for all further testing.
        MultiFab x(ba, dm, Ncomp, Nghost);
        x.setVal(*val);

        Real points = ba.numPts();

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Initial launch to remove any unknown costs in HtoD setup. 

        {
            BL_PROFILE("Initial");
            *val = 0.42;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
            }

            amrex::Print() << "Initial sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }



// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs

        {
            BL_PROFILE("Standard");
            *val = 1.0;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
            }

            amrex::Print() << "No Graph sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create one graph per MFIter iteration and execute them.

        {
            BL_PROFILE("cudaGraph-iter");
            *val = 2.0;

            BL_PROFILE_VAR("cudaGraph-iter-create", cgc);
            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec[x.local_size()];

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Cuda::Device::cudaStream()));

                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("cudaGraph-iter-instantiate", cgi);

            for (int i = 0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("cudaGraph-iter-launch", cgl);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[mfi.LocalIndex()], amrex::Cuda::Device::cudaStream())); 
            }

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Graph-per-iter sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create one graph per CUDA stream and execute them.

        {
            BL_PROFILE("cudaGraph-stream");
            *val = 3.0;

            BL_PROFILE_VAR("cudaGraph-stream-create", cgc);

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

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[i])));
                    }  
                }
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("cudaGraph-stream-instantiate", cgi);

            for (int i = 0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("cudaGraph-stream-launch", cgl);

            for (int i = 0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
            {
                amrex::Gpu::Device::setStreamIndex(i); 
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[i], amrex::Cuda::Device::cudaStream())); 
            }

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();
            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Graph-per-stream sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("cudaGraph-linked-iter");
            *val = 4.0;

            BL_PROFILE_VAR("cudaGraph-linked-iter-create", cgc);

            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Cuda::Device::cudaStream()));

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("cudaGraph-linked-iter-instantiate", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("cudaGraph-linked-iter-launch", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Full-graph-iter sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("cudaGraph-linked-stream");
            *val = 5.0;

            BL_PROFILE_VAR("cudaGraph-linked-stream-create", cgc);

            cudaGraph_t     graph[amrex::Gpu::Device::numCudaStreams()];
            cudaGraphExec_t graphExec;

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

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Cuda::Device::cudaStream(), &(graph[i])));
                    }
                }
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<amrex::Gpu::Device::numCudaStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("cudaGraph-linked-stream-instantiate", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("cudaGraph-linked-stream-launch", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Linked-graph-stream sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;

            BL_PROFILE_VAR("cudaGraph-linked-stream-relaunch", cgrl);

            *val = 10.0;

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Cuda::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();

            BL_PROFILE_VAR_STOP(cgrl);


            amrex::Print() << "Rerun with different val = " << x.sum() << ". Expected = " << points*(*val) << std::endl;


        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
