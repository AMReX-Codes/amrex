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
    amrex::Initialize(argc, argv);
    amrex::Gpu::GraphSafeGuard gpu_gsg(true);
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

        amrex::Print() << "Number of boxes = " << x.size() << std::endl << std::endl;

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
            BL_PROFILE("IterPerGraph");
            *val = 2.0;

            BL_PROFILE_VAR("CREATE: IterPerGraph", cgc);
            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec[x.local_size()];

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                           cudaStreamCaptureModeGlobal));
#endif

                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("INSTANTIATE: IterPerGraph", cgi);

            for (int i = 0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: IterPerGraph", cgl);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[mfi.LocalIndex()], amrex::Gpu::Device::cudaStream())); 
            }

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Graph-per-iter sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create one graph per CUDA stream and execute them.

        {
            BL_PROFILE("StreamPerGraph");
            *val = 3.0;

            BL_PROFILE_VAR("CREATE: StreamPerGraph", cgc);

            cudaGraph_t     graph[amrex::Gpu::numGpuStreams()];
            cudaGraphExec_t graphExec[amrex::Gpu::numGpuStreams()];

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                       AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                       AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                                  cudaStreamCaptureModeGlobal));
#endif
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
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[i])));
                    }  
                }
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("INSTANTIATE: StreamPerGraph", cgi);

            for (int i = 0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: StreamPerGraph", cgl);

            for (int i = 0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                amrex::Gpu::Device::setStreamIndex(i); 
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[i], amrex::Gpu::Device::cudaStream())); 
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
            BL_PROFILE("IterGraph");
            *val = 4.0;

            BL_PROFILE_VAR("CREATE: IterGraph", cgc);

            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                           cudaStreamCaptureModeGlobal));
#endif

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
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
            BL_PROFILE_VAR("INSTANTIATE: IterGraph", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: IterGraph", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            AMREX_GPU_SAFE_CALL(cudaGraphDestroy(graphFull));
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Full-graph-iter sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("StreamGraph");
            *val = 5.0;

            BL_PROFILE_VAR("CREATE: StreamGraph", cgc);

            cudaGraph_t     graph[amrex::Gpu::numGpuStreams()];
            cudaGraphExec_t graphExec;
            cudaEvent_t memcpy_event = {0};
            cudaEventCreate(&memcpy_event);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                        AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                        AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                                   cudaStreamCaptureModeGlobal));
#endif
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
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[i])));
                    }
                }
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("INSTANTIATE: StreamGraph", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: StreamGraph", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            AMREX_GPU_SAFE_CALL(cudaGraphDestroy(graphFull));

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Linked-graph-stream sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;

        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("Indep - Streamless");
            *val = 6.0;

            BL_PROFILE_VAR("CREATE: Indep - Streamless", cgc);

            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                amrex::Gpu::Device::setStreamIndex(0);
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                           cudaStreamCaptureModeGlobal));
#endif

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
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
            BL_PROFILE_VAR("INSTANTIATE: Indep - Streamless", cgi);
/*
            constexpr int log_size = 1028;
            char graph_log[log_size];
            graph_log[0]='\0';

            cudaGraphInstantiate(&graphExec, graphFull, NULL, &(graph_log[0]), log_size);

            amrex::Print() << graph_log << std::endl;
            AMREX_GPU_ERROR_CHECK();
*/
            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: Indep - Streamless", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            AMREX_GPU_SAFE_CALL(cudaGraphDestroy(graphFull));
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();


            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Full-graph-independent-streamless sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("Depend - Streamless");
            *val = 7.0;

            BL_PROFILE_VAR("CREATE: Depend - Streamless", cgc);

            cudaGraph_t     graph;
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                amrex::Gpu::Device::setStreamIndex(0);
                if (mfi.LocalIndex() == 0)
                {
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
                    AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream()));
#else
                    AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(),
                                                               cudaStreamCaptureModeGlobal));
#endif
                } 

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    amrex::Gpu::Device::setStreamIndex(0); 
                    AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &graph));
                }
            }

            BL_PROFILE_VAR_STOP(cgc);
            BL_PROFILE_VAR("INSTANTIATE: Depend - Streamless", cgi);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi);
            BL_PROFILE_VAR("LAUNCH: Depend - Streamless", cgl);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl);

            amrex::Print() << "Linked-graph-dependent-streamless sum = " << x.sum() << ". Expected = " << points*(*val) << std::endl;

        }


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
