#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// TestLoop
// Written as a seperate function for easy changes/testing.

__global__ void fillerKernel() { }

void TestLoopFunc(long* a, int num)
{
    amrex::ParallelFor(num,
    [=] AMREX_GPU_HOST_DEVICE (int i)
    {
        *(a+i) += 1;
    });
}

void Check(AsyncArray<long> &arr, Vector<long> &vec, int value)
{
    arr.copyToHost(vec.data(), vec.size());
    for (int i=0; i<vec.size(); ++i)
    {
        if (vec[i] != (value))
        {
            amrex::Print() << "vec[" << i << "] = " << vec[i]
                           << " != " << value << std::endl;
        }
    }
}

void WriteTimers(std::ofstream& file, long Nnodes,
                 double loop_1, double loop_2, 
                 double buildg_1, double buildg_2, 
                 double instg_1, double instg_2,
                 double rung_1, double rung_2,
                 double graph_1, double graph_2)
{
   file << Nnodes  << " " << loop_1 << " " << buildg_1 << " " << instg_1 << " "
        << rung_1  << " " << graph_1 << " " << loop_2 << " " << buildg_2 << " "
        << instg_2 << " " << rung_2 << " " << graph_2 << std::endl;
}

void InitGraph(int Nnodes)
{
    BL_PROFILE("InitGraph");

    int streams = 16;

    amrex::Print() << "Instantiating a " << Nnodes << " node graph." << std::endl;

    cudaGraphExec_t graphExec;

    for (int n=0; n<(Nnodes); ++n)
    {
        Gpu::Device::startGraphRecording((n == 0), NULL, NULL, 0);

        // ..................
        Gpu::Device::setStreamIndex(n%streams);
        fillerKernel<<<1, 1, 0, Gpu::gpuStream()>>>();
        // ..................

        graphExec = Gpu::Device::stopGraphRecording((n == (Nnodes-1)));
    }

    AMREX_GPU_SAFE_CALL(cudaGraphExecDestroy(graphExec));
}

double Loop(long* h_ptr, long* d_ptr, int Nnodes, int streams, int warp, std::string label)
{
    BL_PROFILE(label);
    double timer = amrex::second(); 

    Gpu::Device::setStreamIndex(0);
    cudaMemcpy(d_ptr, h_ptr, streams*warp, cudaMemcpyHostToDevice);
    for (int n=0; n<(Nnodes*streams); ++n)
    {
        Gpu::Device::setStreamIndex(n%streams);
        size_t offset = (n%streams)*warp;
        TestLoopFunc((d_ptr+offset), warp);
    }
    Gpu::Device::synchronize();
    Gpu::Device::resetStreamIndex();

    timer = amrex::second() - timer; 

    return timer;
}

void Graph(long* h_ptr, long* d_ptr, int Nnodes, int streams, int warp, 
           double& build_time, double& inst_time, double& run_time, double& total_time,
           std::string label)
{
    BL_PROFILE(label);
    total_time = amrex::second();

    BL_PROFILE_VAR("CREATE: " + label, cgc);
    build_time = amrex::second();

    cudaGraph_t     graph;
    cudaGraphExec_t graphExec;

    for (int n=0; n<(Nnodes*streams); ++n)
    {
        if (n == 0)
        {
            Gpu::Device::setStreamIndex(0);
            cudaStream_t graph_stream = Gpu::gpuStream();
            cudaEvent_t memcpy_event = {0};
            AMREX_GPU_SAFE_CALL(cudaEventCreateWithFlags(&memcpy_event, cudaEventDisableTiming));

            AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(graph_stream, cudaStreamCaptureModeGlobal));

            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(d_ptr, h_ptr, warp*streams, cudaMemcpyHostToDevice, graph_stream));
            AMREX_GPU_SAFE_CALL(cudaEventRecord(memcpy_event, graph_stream));

            for (int i=1; i<streams; ++i)
            {
                Gpu::Device::setStreamIndex(i);
                AMREX_GPU_SAFE_CALL(cudaStreamWaitEvent(Gpu::gpuStream(), memcpy_event, 0));
            }

            AMREX_GPU_SAFE_CALL(cudaEventDestroy(memcpy_event));
        }


        // ..................
        Gpu::Device::setStreamIndex(n%streams);
        size_t offset = (n%streams)*warp;
        TestLoopFunc((d_ptr+offset), warp);
        // ..................

        if (n == (Nnodes*streams-1))
        { 
            Gpu::Device::setStreamIndex(0);
            cudaStream_t graph_stream = Gpu::gpuStream();
            cudaEvent_t rejoin_event = {0};
            AMREX_GPU_SAFE_CALL(cudaEventCreateWithFlags(&rejoin_event, cudaEventDisableTiming));

            for (int i=1; i<streams; ++i)
            {
                Gpu::Device::setStreamIndex(i);
                cudaEventRecord(rejoin_event, Gpu::gpuStream());
                cudaStreamWaitEvent(graph_stream, rejoin_event, 0);
            }

            Gpu::Device::resetStreamIndex();

            AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(graph_stream, &graph));
            AMREX_GPU_SAFE_CALL(cudaEventDestroy(rejoin_event));
        }
    }

    build_time = amrex::second() - build_time;
    BL_PROFILE_VAR_STOP(cgc);
    BL_PROFILE_VAR("INSTANTIATE: " + label, cgi);
    inst_time = amrex::second();

    AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));

    inst_time = amrex::second() - inst_time;
    BL_PROFILE_VAR_STOP(cgi);
    BL_PROFILE_VAR("LAUNCH: " + label, cgl);
    run_time = amrex::second();

    amrex::Gpu::Device::setStreamIndex(0); 
    AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 

    amrex::Gpu::Device::synchronize();
    amrex::Gpu::Device::resetStreamIndex();

    run_time = amrex::second() - run_time;
    BL_PROFILE_VAR_STOP(cgl);
    BL_PROFILE_VAR("DESTROY: " + label, cgd);

    AMREX_GPU_SAFE_CALL(cudaGraphExecDestroy(graphExec));
    AMREX_GPU_SAFE_CALL(cudaGraphDestroy(graph));

    BL_PROFILE_VAR_STOP(cgl);

    total_time = amrex::second() - total_time;
}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::Gpu::GraphSafeGuard gsf(true);

        int begin_nodes, end_nodes, streams, factor, warp, init, tofile;
        std::string filename;
        {
            ParmParse pp;

            begin_nodes = 1;
            pp.query("begin_nnodes", begin_nodes);

            end_nodes = 10;
            pp.query("end_nnodes", end_nodes);

            warp = Gpu::Device::warp_size;
            pp.query("threads_per_stream", warp);

            streams = Gpu::Device::numGpuStreams();
            pp.query("nstreams", streams);

            factor = 10;
            pp.query("init_factor", factor);

            filename = "timers.txt";
            pp.query("output", filename);

            init = 0;
            pp.query("init", init);

            tofile = 1;
            pp.query("write_to_file", tofile);
        }

        // Setup data to work on:
        int size = warp*streams;
        Vector<long> vec(size);
        for (int i=0; i<size; ++i)
        {
            vec[i] = 0;
        }
        Gpu::AsyncArray<long> arr(vec.data(), size);
        long* ptr = arr.data();

        if (init)
        {
            InitGraph(end_nodes*16*10);
        }

        amrex::Print() << "Init Graph = " << init << std::endl;
        amrex::Print() << "Number of streams = " << streams << std::endl;
        amrex::Print() << "Threads per stream = " << warp << std::endl;
        amrex::Print() << "Stating number of nodes per stream = " << begin_nodes << std::endl; 
        amrex::Print() << "Ending number of nodes per stream = " << end_nodes << std::endl << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        std::ofstream file;
        if (tofile)
        {
            file.open(filename);
            file.setf(std::ios_base::scientific, std::ios_base::floatfield);
            file << "Nodes / Loop1 / Build1 / Instantiate1 / Run1 / Total1 "
                         << "Loop2 / Build2 / Instantiate2 / Run2 / Total2" << std::endl;
        }

        double loop_1 = 0.0,   loop_2 = 0.0, 
               buildg_1 = 0.0, buildg_2 = 0.0, 
               instg_1 = 0.0,  instg_2 = 0.0, 
               rung_1 = 0.0,   rung_2 = 0.0, 
               graph_1 = 0.0,  graph_2 = 0.0;

        int kidx = 0;

        for (int Nnodes=begin_nodes; Nnodes<=end_nodes; ++Nnodes)
        {
            amrex::Print() << "Testing " << Nnodes << " per stream." << std::endl;

            loop_1 = Loop(vec.data(), ptr, Nnodes, streams, warp, "Loop 1");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            loop_2 = Loop(vec.data(), ptr, Nnodes, streams, warp, "Loop 2");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            Graph(vec.data(), ptr, Nnodes, streams, warp, 
                  buildg_1, instg_1, rung_1, graph_1,
                  "Graph 1");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            loop_2 = Loop(vec.data(), ptr, Nnodes, streams, warp, "Loop 3");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            Graph(vec.data(), ptr, Nnodes, streams, warp, 
                  buildg_2, instg_2, rung_2, graph_2,
                  "Graph 2");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            loop_2 = Loop(vec.data(), ptr, Nnodes, streams, warp, "Loop 4");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            Graph(vec.data(), ptr, Nnodes, streams, warp, 
                  buildg_2, instg_2, rung_2, graph_2,
                  "Graph 3");
            kidx += Nnodes;
            Check(arr, vec, kidx);

            arr.clear();

            if (tofile)
            {
                WriteTimers(file, (Nnodes*streams),
                            loop_1, loop_2, buildg_1, buildg_2, instg_1,
                            instg_2, rung_1, rung_2, graph_1, graph_2);
            }
        }

        if (tofile)
        {
           file.close();
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    } // amrex scope

    amrex::Finalize();
}
