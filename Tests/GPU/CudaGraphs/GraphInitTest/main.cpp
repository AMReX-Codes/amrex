#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// TestLoop
// Written as a seperate function for easy changes/testing.
void TestLoopFunc(long* a, int num)
{
    amrex::ParallelFor(num,
    [=] AMREX_GPU_HOST_DEVICE (int i)
    {
        *(a+i) += 1;
    });
}

void Check(AsyncArray<long> &arr, Vector<long> &vec, int n, int Nnodes)
{
    arr.copyToHost(vec.data(), vec.size());

    for (int i=0; i<vec.size(); ++i)
    {
        if (vec[i] != (n*Nnodes))
        {
            amrex::Print() << "vec[" << i << "] = " << vec[i]
                           << " != " << n*Nnodes << std::endl;
        }
    }
}


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int Nnodes, Nstreams;
        {
            ParmParse pp;

            Nnodes = 1000;
            pp.query("nnodes", Nnodes);

            Nstreams = Gpu::Device::numGpuStreams();
            pp.query("nstreams", Nstreams);
        }

        // Data to "work" on.
        // Assume one warp per launch,
        //   separate workspaces for each stream.
        int warp = Gpu::Device::warp_size;
        int streams = Gpu::Device::numGpuStreams();
        int size = warp*streams;

        Vector<long> vec(size);
        for (int i=0; i<size; ++i)
        {
            vec[i] = 0;
        }

        Gpu::AsyncArray<long> arr(vec.data(), size);
        long* ptr = arr.data();

        amrex::Print() << "Number of points = "   << size << std::endl;
        amrex::Print() << "Number of launches = " << Nnodes << std::endl << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Initial launch to remove any unknown costs in HtoD setup. 

        {
            amrex::Gpu::GraphSafeGuard gsg(false);

            BL_PROFILE("Initial");

            for (int n=0; n<(Nnodes*streams); ++n)
            {
                Gpu::Device::setStreamIndex(n%streams);
                size_t offset = (n%streams)*warp;
                TestLoopFunc((ptr+offset), warp);
            }
            Gpu::Device::synchronize();
            Gpu::Device::resetStreamIndex();
        }

        Check(arr, vec, 1, Nnodes);

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs

        {
            amrex::Gpu::GraphSafeGuard gsg(false);

            BL_PROFILE("Loop w/ Streams");

            for (int n=0; n<(Nnodes*streams); ++n)
            {
                Gpu::Device::setStreamIndex(n%streams);
                size_t offset = (n%streams)*warp;
                TestLoopFunc((ptr+offset), warp);
            }
            Gpu::Device::synchronize();
            Gpu::Device::resetStreamIndex();
        }

        Check(arr, vec, 2, Nnodes);


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("StreamGraph");

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
                        AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));
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

    }

    amrex::Finalize();
}
