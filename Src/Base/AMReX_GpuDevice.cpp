
#include <iostream>
#include <map>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <AMReX_GpuDevice.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_GpuLaunch.H>

#if defined(AMREX_USE_CUDA)
#include <cuda_profiler_api.h>
#if defined(AMREX_PROFILING) || defined (AMREX_TINY_PROFILING)
#include "nvToolsExt.h"
#endif
#endif

#ifdef AMREX_USE_ACC
extern "C" {
    void amrex_initialize_acc (int);
    void amrex_finalize_acc ();
}
#endif

namespace amrex {
namespace Gpu {

int Device::device_id = 0;
int Device::num_devices_used = 0;
int Device::verbose = 0;

#ifdef AMREX_USE_GPU
constexpr int Device::max_gpu_streams;
dim3 Device::numThreadsMin      = dim3(1, 1, 1);
dim3 Device::numThreadsOverride = dim3(0, 0, 0);
dim3 Device::numBlocksOverride  = dim3(0, 0, 0);
int  Device::max_blocks_per_launch = 640;

std::array<gpuStream_t,Device::max_gpu_streams> Device::gpu_streams;
gpuStream_t                                     Device::gpu_stream;
gpuDeviceProp_t                                 Device::device_prop;

constexpr int                                   Device::warp_size;

namespace {

    AMREX_GPU_GLOBAL void emptyKernel() {}

    void InitializeGraph(int graph_size)
    {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )

        BL_PROFILE("InitGraph");

        int streams = Gpu::Device::numGpuStreams();
        cudaGraphExec_t graphExec;
        for (int n=0; n<(graph_size); ++n)
        {
            Gpu::Device::startGraphRecording((n == 0), NULL, NULL, 0);

            // ..................
            Gpu::Device::setStreamIndex(n%streams);
            emptyKernel<<<1, 1, 0, Gpu::gpuStream()>>>();
            // ..................

            graphExec = Gpu::Device::stopGraphRecording((n == (graph_size-1)));
        }
        AMREX_CUDA_SAFE_CALL(cudaGraphExecDestroy(graphExec));
#endif
    }
}

#endif

void
Device::Initialize ()
{
#ifdef AMREX_USE_HIP

    ParmParse pp("device");

    pp.query("v", verbose);
    pp.query("verbose", verbose);

    if (amrex::Verbose()) {
        amrex::Print() << "Initializing HIP...\n";
    }

    int gpu_device_count;
    AMREX_HIP_SAFE_CALL(hipGetDeviceCount(&gpu_device_count));

    if (gpu_device_count <= 0) {
        amrex::Abort("No GPU device found");
    }

    // Now, assign ranks to GPUs. If we only have one GPU,
    // or only one MPI rank, this is easy. Otherwise, we
    // need to do a little more work.

    if (ParallelDescriptor::NProcs() == 1) {
        device_id = 0;
    }
    else if (gpu_device_count == 1) {
        device_id = 0;
    }
    else {
        amrex::Abort("USE_HIP and USE_MPI not supported yet");
    }

    AMREX_HIP_SAFE_CALL(hipSetDevice(device_id));
    AMREX_HIP_SAFE_CALL(hipSetDeviceFlags(hipDeviceMapHost));

    initialize_gpu();

    if (amrex::Verbose()) {
#ifdef AMREX_USE_MPI
        amrex::Print() << "HIP initialized with 1 GPU per MPI rank\n";
#else
        amrex::Print() << "HIP initialized with 1 GPU\n";
#endif
    }

#elif defined(AMREX_USE_CUDA)

#if defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING)
    // Wrap cuda init to identify it appropriately in nvvp.
    // Note: first substantial cuda call may cause a lengthy
    // cuda API and cuda driver API initialization that will
    // be captured by the profiler. It a necessary, system
    // dependent step that is unavoidable.
    nvtxRangeId_t nvtx_init;
    const char* pname = "initialize_device";
    nvtx_init = nvtxRangeStartA(pname);
#endif

    ParmParse pp("device");

    pp.query("v", verbose);
    pp.query("verbose", verbose);

    if (amrex::Verbose()) {
        amrex::Print() << "Initializing CUDA...\n";
    }

    // XL CUDA Fortran support needs to be initialized
    // before any CUDA API calls.

#if (defined(__ibmxl__) && !defined(BL_NO_FORT))
    __xlcuf_init();
#endif

    // Count the number of CUDA visible devices.

    int cuda_device_count;
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceCount(&cuda_device_count));

    if (cuda_device_count <= 0) {
        amrex::Abort("No CUDA device found");
    }

    // Now, assign ranks to GPUs. If we only have one GPU,
    // or only one MPI rank, this is easy. Otherwise, we
    // need to do a little more work.

    if (ParallelDescriptor::NProcs() == 1) {
        device_id = 0;
    }
    else if (cuda_device_count == 1) {
        device_id = 0;
    }
    else {

        // ifdef the following against MPI so it compiles, but note
        // that we can only get here if using more than one processor,
        // which requires MPI.

#ifdef BL_USE_MPI

        // Create a communicator out of only the ranks sharing GPUs.
        // The default assumption is that this is all the ranks on the
        // same node, and to get that we'll use the MPI-3.0 split that
        // looks for shared memory communicators (and we'll error out
        // if that standard is unsupported).

#if MPI_VERSION < 3
        amrex::Abort("When using CUDA with MPI, if multiple devices are visible to each rank, MPI-3.0 must be supported.");
#endif

        // However, it's possible that the ranks sharing GPUs will be
        // confined to a single socket rather than a full node. Indeed,
        // this is often the optimal configuration; for example, on Summit,
        // a good configuration using jsrun is one resource set per
        // socket (two per node), with three GPUs per resource set.
        // To deal with this where we can, we'll take advantage of OpenMPI's
        // specialized split by socket. However, we only want to do this
        // if in fact our resource set is confined to the socket.
        // To make this determination we need to have system information,
        // which is provided by the build system for the systems
        // we know about. The simple heuristic we'll use to determine
        // this is if the number of visible devices is smaller than
        // the known number of GPUs per socket.

#if (!defined(AMREX_GPUS_PER_SOCKET) && !defined(AMREX_GPUS_PER_NODE))
        amrex::Warning("Multiple GPUs are visible to each MPI rank, but the number of GPUs per socket or node has not been provided.\n"
                       "This may lead to incorrect or suboptimal rank-to-GPU mapping.");
#endif

        MPI_Comm local_comm;

        int split_type;

#if (defined(OPEN_MPI) && defined(AMREX_GPUS_PER_SOCKET))
        if (cuda_device_count <= AMREX_GPUS_PER_SOCKET)
            split_type = OMPI_COMM_TYPE_SOCKET;
        else
            split_type = OMPI_COMM_TYPE_NODE;
#else
        split_type = MPI_COMM_TYPE_SHARED;
#endif

        // We have no preference on how ranks get ordered within this communicator.
        int key = 0;

        MPI_Comm_split_type(ParallelDescriptor::Communicator(), split_type, key, MPI_INFO_NULL, &local_comm);

        // Get rank within the local communicator, and number of ranks.
        int n_procs;
        MPI_Comm_size(local_comm, &n_procs);

        int my_rank;
        MPI_Comm_rank(local_comm, &my_rank);

        // Free the local communicator.
        MPI_Comm_free(&local_comm);

        // For each rank that shares a GPU, use round-robin assignment
        // to assign MPI ranks to GPUs. We will arbitrarily assign
        // ranks to GPUs, assuming that socket awareness has already
        // been handled.

        device_id = my_rank % cuda_device_count;

        // If we detect more ranks than visible GPUs, warn the user
        // that this will fail in the case where the devices are
        // set to exclusive process mode and MPS is not enabled.

        if (n_procs > cuda_device_count) {
            amrex::Print() << "Mapping more than one rank per GPU. This will fail if the GPUs are in exclusive process mode\n"
                           << "and MPS is not enabled. In that case you will see an error such as all CUDA-capable devices are\n"
                           << "busy. To resolve that issue, set the GPUs to the default compute mode, or enable MPS. If you are\n"
                           << "on a cluster, please consult the system user guide for how to launch your job in this configuration.\n";
        }

#endif

    }

    AMREX_CUDA_SAFE_CALL(cudaSetDevice(device_id));
    AMREX_CUDA_SAFE_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));

#ifdef AMREX_USE_ACC
    amrex_initialize_acc(device_id);
#endif

    initialize_gpu();

    // Count up the total number of devices used by
    // all MPI ranks. Since we have to consider the
    // case of multiple ranks per GPU, we cannot simply
    // set it to the number of MPI ranks. A reliable way
    // to do this instead is to collect the UUID of each
    // GPU used by every rank, perform a gather, and then
    // count the number of unique UUIDs in the result.

    // Note: the field we need from the CUDA device properties
    // is only available starting from CUDA 10.0, so we will
    // leave num_devices_used as 0 for older CUDA toolkits.

#if AMREX_NVCC_MAJOR_VERSION >= 10
    size_t uuid_length = 16;
    size_t recv_sz = uuid_length * ParallelDescriptor::NProcs();
    const char* sendbuf = &device_prop.uuid.bytes[0];
    char* recvbuf = new char[recv_sz];

    ParallelDescriptor::Gather<char,char>(sendbuf, uuid_length,
                                          recvbuf, uuid_length,
                                          ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
        std::unordered_set<std::string> uuids;
        for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
            std::string uuid(&recvbuf[16 * i], 16);
            if (uuids.find(uuid) == uuids.end()) {
                uuids.insert(uuid);
            }
        }
        num_devices_used = uuids.size();
    }
    ParallelDescriptor::Bcast<int>(&num_devices_used, 1);

    delete[] recvbuf;
#endif

#if defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING)
    nvtxRangeEnd(nvtx_init);
#endif

    if (amrex::Verbose()) {
#if defined(AMREX_USE_MPI) && (AMREX_NVCC_MAJOR_VERSION >= 10)
        amrex::Print() << "CUDA initialized with 1 GPU per MPI rank; "
                       << num_devices_used << " GPU(s) used in total\n";
#else
        amrex::Print() << "CUDA initialized with 1 GPU\n";
#endif
    }

    cudaProfilerStart();

#endif

}

void
Device::Finalize ()
{
#if defined(AMREX_USE_HIP)

    for (int i = 0; i < max_gpu_streams; ++i)
    {
        AMREX_HIP_SAFE_CALL(hipStreamDestroy(gpu_streams[i]));
    }

    AMREX_HIP_SAFE_CALL(hipDeviceReset());

#elif defined(AMREX_USE_CUDA)

    cudaProfilerStop();

    for (int i = 0; i < max_gpu_streams; ++i)
    {
        AMREX_CUDA_SAFE_CALL(cudaStreamDestroy(gpu_streams[i]));
    }

#ifdef AMREX_USE_ACC
    amrex_finalize_acc();
#endif

    AMREX_CUDA_SAFE_CALL(cudaDeviceReset());
#endif
}

void
Device::initialize_gpu ()
{
#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_HIP

    AMREX_HIP_SAFE_CALL(hipGetDeviceProperties(&device_prop, device_id));

    // check compute capability

    if (sizeof(Real) == 8) {
        AMREX_HIP_SAFE_CALL(hipDeviceSetSharedMemConfig(hipSharedMemBankSizeEightByte));
    } else if (sizeof(Real) == 4) {
        AMREX_HIP_SAFE_CALL(hipDeviceSetSharedMemConfig(hipSharedMemBankSizeFourByte));
    }

    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_HIP_SAFE_CALL(hipStreamCreate(&gpu_streams[i]));
    }

#else
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(device_prop.major >= 6, "Compute capability must be >= 6");

    // Prefer L1 cache to shared memory (this has no effect on GPUs with a fixed L1 cache size).
    AMREX_CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    if (sizeof(Real) == 8) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    } else if (sizeof(Real) == 4) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
    }

    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_CUDA_SAFE_CALL(cudaStreamCreate(&gpu_streams[i]));
    }

#endif

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warp_size == device_prop.warpSize, "Incorrect warp size");

    gpu_stream = 0;

    ParmParse pp("device");

    int nx = 0;
    int ny = 0;
    int nz = 0;

    pp.query("numThreads.x", nx);
    pp.query("numThreads.y", ny);
    pp.query("numThreads.z", nz);

    numThreadsOverride.x = (int) nx;
    numThreadsOverride.y = (int) ny;
    numThreadsOverride.z = (int) nz;

    nx = 0;
    ny = 0;
    nz = 0;

    pp.query("numBlocks.x", nx);
    pp.query("numBlocks.y", ny);
    pp.query("numBlocks.z", nz);

    numBlocksOverride.x = (int) nx;
    numBlocksOverride.y = (int) ny;
    numBlocksOverride.z = (int) nz;

    // Graph initialization
    int graph_init = 0;
    int graph_size = 10000;
    pp.query("graph_init", graph_init);
    pp.query("graph_init_nodes", graph_size);

    if (graph_init)
    {
        GraphSafeGuard gsg(true);
        InitializeGraph(graph_size);
    }

    max_blocks_per_launch = numMultiProcessors() * maxThreadsPerMultiProcessor() / AMREX_GPU_MAX_THREADS;

#endif
}

int
Device::deviceId () noexcept
{
    return device_id;
}

int
Device::numDevicesUsed () noexcept
{
    return num_devices_used;
}

void
Device::setStreamIndex (const int idx) noexcept
{
#ifdef AMREX_USE_GPU
    if (idx < 0) {
        gpu_stream = 0;
    } else {
        gpu_stream = gpu_streams[idx % max_gpu_streams];
    }
#endif
}

#ifdef AMREX_USE_GPU
gpuStream_t
Device::resetStream () noexcept
{
    gpuStream_t r = gpu_stream;
    gpu_stream = 0;
    return r;
}

gpuStream_t
Device::setStream (gpuStream_t s) noexcept
{
    gpuStream_t r = gpu_stream;
    gpu_stream = s;
    return r;
}
#endif

void
Device::synchronize () noexcept
{
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipDeviceSynchronize());,
                       AMREX_CUDA_SAFE_CALL(cudaDeviceSynchronize()); )
}

void
Device::streamSynchronize () noexcept
{
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipStreamSynchronize(gpu_stream));,
                       AMREX_CUDA_SAFE_CALL(cudaStreamSynchronize(gpu_stream)); )
}


#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )

void
Device::startGraphRecording(bool first_iter, void* h_ptr, void* d_ptr, size_t sz)
{
    if ((first_iter) && inLaunchRegion() && inGraphRegion())
    {
        // Uses passed information to do initial async memcpy in graph and 
        //    links dependency to all streams using cudaEvents.

        setStreamIndex(0);
        cudaStream_t graph_stream = gpuStream();
        cudaEvent_t memcpy_event = {0};
        AMREX_CUDA_SAFE_CALL( cudaEventCreate(&memcpy_event, cudaEventDisableTiming) );

#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
        AMREX_CUDA_SAFE_CALL(cudaStreamBeginCapture(graph_stream));
#else  
        AMREX_CUDA_SAFE_CALL(cudaStreamBeginCapture(graph_stream, cudaStreamCaptureModeGlobal));
#endif

        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(d_ptr, h_ptr, sz, cudaMemcpyHostToDevice, graph_stream));
        AMREX_CUDA_SAFE_CALL(cudaEventRecord(memcpy_event, graph_stream));

        // Note: Main graph stream fixed at 0, so i starts at 1.
        //       Will need more complex logic if this changes.
        for (int i=1; i<numGpuStreams(); ++i)
        {
            setStreamIndex(i);
            AMREX_CUDA_SAFE_CALL(cudaStreamWaitEvent(gpuStream(), memcpy_event, 0));
        }
        setStreamIndex(0);

        AMREX_CUDA_SAFE_CALL( cudaEventDestroy(memcpy_event) );
    }
}

cudaGraphExec_t
Device::stopGraphRecording(bool last_iter)
{
    cudaGraphExec_t graphExec;

    if (last_iter && inLaunchRegion() && inGraphRegion())
    {
        // Uses cudaEvents to rejoin the streams, making a single graph.
        setStreamIndex(0);
        cudaStream_t graph_stream = gpuStream();
        cudaEvent_t rejoin_event = {0};
        AMREX_CUDA_SAFE_CALL( cudaEventCreate(&rejoin_event, cudaEventDisableTiming) );

        // Note: Main graph stream fixed at 0, so i starts at 1.
        //       Will need more complex logic if this changes.
        for (int i=1; i<Gpu::Device::numGpuStreams(); ++i)
        {
            Gpu::Device::setStreamIndex(i);
            cudaEventRecord(rejoin_event, gpuStream());
            cudaStreamWaitEvent(graph_stream, rejoin_event, 0);
        }
        Gpu::Device::setStreamIndex(0);

        cudaGraph_t graph;
        AMREX_CUDA_SAFE_CALL(cudaStreamEndCapture(graph_stream, &graph));
        graphExec = instantiateGraph(graph);

        AMREX_CUDA_SAFE_CALL( cudaGraphDestroy(graph); );
        AMREX_CUDA_SAFE_CALL( cudaEventDestroy(rejoin_event) );
    }

    return graphExec;
}

cudaGraphExec_t
Device::instantiateGraph(cudaGraph_t graph)
{
    cudaGraphExec_t graphExec;

#ifdef AMREX_DEBUG 
//  Implementes cudaGraphInstantiate error logging feature.
//  Upon error, delays abort until message is output. 
    constexpr int log_size = 1028;
    char graph_log[log_size];
    graph_log[0]='\0';

    cudaGraphInstantiate(&graphExec, graph, NULL, &(graph_log[0]), log_size); 

    if (graph_log[0] != '\0')
    {
        amrex::Print() << graph_log << std::endl;
        AMREX_GPU_ERROR_CHECK();
    }
#else

    AMREX_CUDA_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0)); 

#endif

    return graphExec;

}

void
Device::executeGraph(const cudaGraphExec_t &graphExec, bool synch)
{
    if (inLaunchRegion() && inGraphRegion())
    {
        setStreamIndex(0);
        AMREX_CUDA_SAFE_CALL(cudaGraphLaunch(graphExec, cudaStream()));
        if (synch) {
            synchronize();
        }
        resetStreamIndex();
    }
}

#endif

void
Device::mem_advise_set_preferred (void* p, const std::size_t sz, const int device)
{
    // HIP does not support memory advise.
#ifdef AMREX_USE_CUDA
#ifndef AMREX_USE_HIP
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
#endif
    {
        AMREX_CUDA_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device));
    }
#endif

}

void
Device::mem_advise_set_readonly (void* p, const std::size_t sz)
{
    // HIP does not support memory advise.
#ifdef AMREX_USE_CUDA
#ifndef AMREX_USE_HIP
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
#endif
    {
        AMREX_CUDA_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId));
    }
#endif
}

#ifdef AMREX_USE_GPU

void
Device::setNumThreadsMin (int nx, int ny, int nz) noexcept
{
    numThreadsMin.x = nx;
    numThreadsMin.y = ny;
    numThreadsMin.z = nz;
}

void
Device::n_threads_and_blocks (const long N, dim3& numBlocks, dim3& numThreads) noexcept
{
    numThreads = AMREX_GPU_MAX_THREADS;
    numBlocks = std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, 1L); // in case N = 0
}

void
Device::c_comps_threads_and_blocks (const int* lo, const int* hi, const int comps,
                                    dim3& numBlocks, dim3& numThreads) noexcept
{
    c_threads_and_blocks(lo, hi, numBlocks, numThreads);
    numBlocks.x *= static_cast<unsigned>(comps);
}

void
Device::c_threads_and_blocks (const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) noexcept
{
    // Our threading strategy will be to allocate thread blocks
    //preferring the x direction first to guarantee coalesced accesses.
    int tile_size[] = {AMREX_D_DECL(hi[0]-lo[0]+1,hi[1]-lo[1]+1,hi[2]-lo[2]+1)};

#if (AMREX_SPACEDIM == 1)

    numThreads.x = std::min(tile_size[0], AMREX_GPU_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = 1;
    numBlocks.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(tile_size[0]), AMREX_GPU_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(tile_size[1]), AMREX_GPU_MAX_THREADS / numThreads.x   );
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = (tile_size[1] + numThreads.y - 1) / numThreads.y;
    numBlocks.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, std::min(static_cast<unsigned>(tile_size[0]), numThreads.x));
    numThreads.y = std::max(numThreadsMin.y, std::min(static_cast<unsigned>(tile_size[1]), numThreads.y));
    numThreads.z = std::max(numThreadsMin.z, std::min(static_cast<unsigned>(tile_size[2]), numThreads.z));

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = (tile_size[1] + numThreads.y - 1) / numThreads.y;
    numBlocks.z = (tile_size[2] + numThreads.z - 1) / numThreads.z;

#endif

    AMREX_ASSERT(numThreads.x <= device_prop.maxThreadsDim[0]);
    AMREX_ASSERT(numThreads.y <= device_prop.maxThreadsDim[1]);
    AMREX_ASSERT(numThreads.z <= device_prop.maxThreadsDim[2]);
    AMREX_ASSERT(numThreads.x*numThreads.y*numThreads.z <= device_prop.maxThreadsPerBlock);
    AMREX_ASSERT(numThreads.x > 0);
    AMREX_ASSERT(numThreads.y > 0);
    AMREX_ASSERT(numThreads.z > 0);
    AMREX_ASSERT(numBlocks.x > 0);
    AMREX_ASSERT(numBlocks.y > 0);
    AMREX_ASSERT(numBlocks.z > 0);
}

void
Device::grid_stride_threads_and_blocks (dim3& numBlocks, dim3& numThreads) noexcept
{
    int num_SMs = device_prop.multiProcessorCount;

    int SM_mult_factor = 32;

    if (num_SMs > 0) {

        numBlocks.x = 1;
        numBlocks.y = SM_mult_factor;
        numBlocks.z = num_SMs;

    } else {

        // Arbitrarily set this to a somewhat large number.

        numBlocks.x = 1000;
        numBlocks.y = 1;
        numBlocks.z = 1;

    }

#if (AMREX_SPACEDIM == 1)

    numThreads.x = std::min(device_prop.maxThreadsDim[0], AMREX_GPU_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / numThreads.x);
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = std::max(numThreadsMin.z, numThreads.z);

#endif

    // Allow the user to override these at runtime.

    if (numBlocksOverride.x > 0)
        numBlocks.x = numBlocksOverride.x;
    if (numBlocksOverride.y > 0)
        numBlocks.y = numBlocksOverride.y;
    if (numBlocksOverride.z > 0)
        numBlocks.z = numBlocksOverride.z;

    if (numThreadsOverride.x > 0)
        numThreads.x = numThreadsOverride.x;
    if (numThreadsOverride.y > 0)
        numThreads.y = numThreadsOverride.y;
    if (numThreadsOverride.z > 0)
        numThreads.z = numThreadsOverride.z;

}

void
Device::box_threads_and_blocks (const Box& bx, dim3& numBlocks, dim3& numThreads) noexcept
{
    int num_SMs = device_prop.multiProcessorCount;

    int SM_mult_factor = 32;

    if (num_SMs > 0) {

        numBlocks.x = 1;
        numBlocks.y = SM_mult_factor;
        numBlocks.z = num_SMs;

    } else {

        // Arbitrarily set this to a somewhat large number.

        numBlocks.x = 1000;
        numBlocks.y = 1;
        numBlocks.z = 1;

    }

#if (AMREX_SPACEDIM == 1)

    numThreads.x = std::min(device_prop.maxThreadsDim[0], AMREX_GPU_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

    // Limit the number of threads per block to be no larger in each dimension
    // than the sizes of the box in the corresponding dimension.
    numThreads.x = std::min(numThreads.x, static_cast<unsigned>(bx.length(0)));
    
#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / numThreads.x);
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

    // Limit the number of threads per block to be no larger in each dimension
    // than the sizes of the box in the corresponding dimension.
    numThreads.x = std::min(numThreads.x, static_cast<unsigned>(bx.length(0)));
    numThreads.y = std::min(numThreads.y, static_cast<unsigned>(bx.length(1)));
    
#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = std::max(numThreadsMin.z, numThreads.z);

    // Limit the number of threads per block to be no larger in each dimension
    // than the sizes of the box in the corresponding dimension.
    numThreads.x = std::min(numThreads.x, static_cast<unsigned>(bx.length(0)));
    numThreads.y = std::min(numThreads.y, static_cast<unsigned>(bx.length(1)));
    numThreads.z = std::min(numThreads.z, static_cast<unsigned>(bx.length(2)));
    
#endif

    // Allow the user to override these at runtime.

    if (numBlocksOverride.x > 0)
        numBlocks.x = numBlocksOverride.x;
    if (numBlocksOverride.y > 0)
        numBlocks.y = numBlocksOverride.y;
    if (numBlocksOverride.z > 0)
        numBlocks.z = numBlocksOverride.z;

    if (numThreadsOverride.x > 0)
        numThreads.x = numThreadsOverride.x;
    if (numThreadsOverride.y > 0)
        numThreads.y = numThreadsOverride.y;
    if (numThreadsOverride.z > 0)
        numThreads.z = numThreadsOverride.z;

}

#endif

std::size_t
Device::freeMemAvailable ()
{
#ifdef AMREX_USE_GPU
    std::size_t f, t;
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipMemGetInfo(&f,&t));,
                       AMREX_CUDA_SAFE_CALL(cudaMemGetInfo(&f,&t)); )
    return f;
#else
    return 0;
#endif
}

#ifdef AMREX_USE_GPU
namespace {
    static int ncallbacks = 0;
}

void callbackAdded ()
{
    ++ncallbacks;
}

void resetNumCallbacks ()
{
    ncallbacks = 0;
}

int getNumCallbacks ()
{
    return ncallbacks;
}
#endif

}}
