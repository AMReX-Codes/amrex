
#include <iostream>
#include <map>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <exception>
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
#include <openacc.h>

extern "C" {
    void amrex_initialize_acc (int);
    void amrex_finalize_acc ();
    void amrex_set_acc_stream (int);
}
#endif

#ifdef AMREX_USE_DPCPP
namespace {
    auto amrex_sycl_error_handler = [] (sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch (sycl::exception const& ex) {
                amrex::Abort(std::string("Async SYCL exception: ")+ex.what()+"!!!!!");
            }
        }
    };
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
gpuStream_t                                     Device::gpu_default_stream;
Vector<gpuStream_t>                             Device::gpu_stream;
gpuDeviceProp_t                                 Device::device_prop;

constexpr int                                   Device::warp_size;

#ifdef AMREX_USE_DPCPP
std::unique_ptr<sycl::context> Device::sycl_context;
std::unique_ptr<sycl::device>  Device::sycl_device;
#endif

namespace {

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
    AMREX_GPU_GLOBAL void emptyKernel() {}
#endif

    void InitializeGraph(int graph_size)
    {
        amrex::ignore_unused(graph_size);

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

#if defined(AMREX_USE_CUDA) && (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
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
        AMREX_HIP_OR_CUDA_OR_DPCPP
            ( amrex::Print() << "Initializing HIP...\n";,
              amrex::Print() << "Initializing CUDA...\n";,
              amrex::Print() << "Initializing oneAPI...\n"; )
    }

    // XL CUDA Fortran support needs to be initialized
    // before any CUDA API calls.

#if (defined(AMREX_USE_CUDA) && defined(__ibmxl__) && !defined(BL_NO_FORT))
    __xlcuf_init();
#endif

    // Count the number of GPU devices.
    int gpu_device_count = 0;
#ifdef AMREX_USE_DPCPP
    {
        sycl::gpu_selector device_selector;
        sycl::platform platform(device_selector);
        auto const& gpu_devices = platform.get_devices();
        gpu_device_count = gpu_devices.size();
        if (gpu_device_count <= 0) {
            amrex::Abort("No GPU device found");
        } else if (gpu_device_count > 1) {
            amrex::Abort("DPCPP TODO: more than one device not supported yet");
        }
    }
#else
    AMREX_HIP_OR_CUDA(AMREX_HIP_SAFE_CALL (hipGetDeviceCount(&gpu_device_count));,
                      AMREX_CUDA_SAFE_CALL(cudaGetDeviceCount(&gpu_device_count)); );
    if (gpu_device_count <= 0) {
        amrex::Abort("No GPU device found");
    }
#endif

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
        amrex::Abort("When using GPUs with MPI, if multiple devices are visible to each rank, MPI-3.0 must be supported.");
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
        if (gpu_device_count <= AMREX_GPUS_PER_SOCKET)
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

        device_id = my_rank % gpu_device_count;

        // If we detect more ranks than visible GPUs, warn the user
        // that this will fail in the case where the devices are
        // set to exclusive process mode and MPS is not enabled.

        if (n_procs > gpu_device_count) {
            amrex::Print() << "Mapping more than one rank per GPU. This will fail if the GPUs are in exclusive process mode\n"
                           << "and MPS is not enabled. In that case you will see an error such as: 'all CUDA-capable devices are\n"
                           << "busy'. To resolve that issue, set the GPUs to the default compute mode, or enable MPS. If you are\n"
                           << "on a cluster, please consult the system user guide for how to launch your job in this configuration.\n";
        }

#endif   // BL_USE_MPI

    }

    AMREX_HIP_OR_CUDA(AMREX_HIP_SAFE_CALL (hipSetDevice(device_id));,
                      AMREX_CUDA_SAFE_CALL(cudaSetDevice(device_id)); );
    AMREX_HIP_OR_CUDA(AMREX_HIP_SAFE_CALL (hipSetDeviceFlags(hipDeviceMapHost));,
                      AMREX_CUDA_SAFE_CALL(cudaSetDeviceFlags(cudaDeviceMapHost)); );

#ifdef AMREX_USE_ACC
    amrex_initialize_acc(device_id);
#endif

    initialize_gpu();

#ifdef AMREX_USE_CUDA

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

#if (__CUDACC_VER_MAJOR__ >= 10)
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

#if (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
    nvtxRangeEnd(nvtx_init);
#endif
    cudaProfilerStart();

    if (amrex::Verbose()) {
#if defined(AMREX_USE_MPI) && (__CUDACC_VER_MAJOR__ >= 10)
        if (num_devices_used == ParallelDescriptor::NProcs())
        {
            amrex::Print() << "CUDA initialized with 1 GPU per MPI rank; "
                           << num_devices_used << " GPU(s) used in total\n";
        }
        else 
        {
            amrex::Print() << "CUDA initialized with " << num_devices_used << " GPU(s) and "
                           << ParallelDescriptor::NProcs() << " ranks.\n";
        }
#else  // Should always be using NVCC >= 10 now, so not going to bother with other combinations.
        amrex::Print() << "CUDA initialized with 1 GPU\n"; 
#endif // AMREX_USE_MPI && NVCC >= 10
    }

    cudaProfilerStart();

#elif defined(AMREX_USE_HIP) 
    if (amrex::Verbose()) {
        amrex::Print() << "HIP initialized.\n";
    }
#elif defined(AMREX_USE_DPCPP)
    if (amrex::Verbose()) {
        amrex::Print() << "oneAPI initialized.\n";
    }
#endif

}

void
Device::Finalize ()
{
#ifdef AMREX_USE_CUDA
    cudaProfilerStop();
#endif

    for (int i = 0; i < max_gpu_streams; ++i)
    {
        AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL( hipStreamDestroy(gpu_streams[i]));,
                          AMREX_CUDA_SAFE_CALL(cudaStreamDestroy(gpu_streams[i])); );
    }

#ifdef AMREX_USE_DPCPP
    sycl_context.reset();
    sycl_device.reset();
    for (auto& s : gpu_streams) {
        delete s.queue;
        s.queue = nullptr;
    }
    gpu_stream.clear();
    delete gpu_default_stream.queue;
    gpu_default_stream.queue = nullptr;
#endif

#ifdef AMREX_USE_ACC
    amrex_finalize_acc();
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

    gpu_default_stream = 0;
    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_HIP_SAFE_CALL(hipStreamCreate(&gpu_streams[i]));
    }

#elif defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(device_prop.major >= 6, "Compute capability must be >= 6");

    // Prefer L1 cache to shared memory (this has no effect on GPUs with a fixed L1 cache size).
    AMREX_CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    if (sizeof(Real) == 8) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    } else if (sizeof(Real) == 4) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
    }

    gpu_default_stream = 0;
    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_CUDA_SAFE_CALL(cudaStreamCreate(&gpu_streams[i]));
#ifdef AMREX_USE_ACC
        acc_set_cuda_stream(i, gpu_streams[i]);
#endif
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warp_size == device_prop.warpSize, "Incorrect warp size");

#elif defined(AMREX_USE_DPCPP)
    { // create device, context and queues
        sycl::gpu_selector device_selector;
        sycl_device.reset(new sycl::device(device_selector));
        sycl_context.reset(new sycl::context(*sycl_device, amrex_sycl_error_handler));
        gpu_default_stream.queue = new sycl::queue(*sycl_context, device_selector,
                                         sycl::property_list{sycl::property::queue::in_order{}});
        for (int i = 0; i < max_gpu_streams; ++i) {
            gpu_streams[i].queue = new sycl::queue(*sycl_context, device_selector,
                                         sycl::property_list{sycl::property::queue::in_order{}});
        }
    }

    { // device property
        auto const& d = *sycl_device;
        device_prop.name = d.get_info<sycl::info::device::name>();
        device_prop.totalGlobalMem = d.get_info<sycl::info::device::global_mem_size>();
        device_prop.sharedMemPerBlock = d.get_info<sycl::info::device::local_mem_size>();
        device_prop.multiProcessorCount = d.get_info<sycl::info::device::max_compute_units>();
        device_prop.maxThreadsPerMultiProcessor = -1; // xxxxx DPCPP todo: d.get_info<sycl::info::device::max_work_items_per_compute_unit>(); // unknown
        device_prop.maxThreadsPerBlock = d.get_info<sycl::info::device::max_work_group_size>();
        auto mtd = d.get_info<sycl::info::device::max_work_item_sizes>();
        device_prop.maxThreadsDim[0] = mtd[0];
        device_prop.maxThreadsDim[1] = mtd[1];
        device_prop.maxThreadsDim[2] = mtd[2];
        device_prop.maxGridSize[0] = -1; // xxxxx DPCPP todo: unknown
        device_prop.maxGridSize[0] = -1; // unknown
        device_prop.maxGridSize[0] = -1; // unknown
        device_prop.warpSize = warp_size;
        auto sgss = d.get_info<sycl::info::device::sub_group_sizes>();
        device_prop.maxMemAllocSize = d.get_info<sycl::info::device::max_mem_alloc_size>();
        device_prop.managedMemory = d.get_info<sycl::info::device::host_unified_memory>();
        device_prop.concurrentManagedAccess = d.get_info<sycl::info::device::usm_shared_allocations>();
        device_prop.maxParameterSize = d.get_info<sycl::info::device::max_parameter_size>();
        {
            amrex::Print() << "Device Properties:\n"
                           << "  name: " << device_prop.name << "\n"
                           << "  totalGlobalMem: " << device_prop.totalGlobalMem << "\n"
                           << "  sharedMemPerBlock: " << device_prop.sharedMemPerBlock << "\n"
                           << "  multiProcessorCount: " << device_prop.multiProcessorCount << "\n"
                           << "  maxThreadsPerBlock: " << device_prop.maxThreadsPerBlock << "\n"
                           << "  maxThreadsDim: (" << device_prop.maxThreadsDim[0] << ", " << device_prop.maxThreadsDim[1] << ", " << device_prop.maxThreadsDim[2] << ")\n"
                           << "  warpSize:";
            for (auto s : sgss) {
                amrex::Print() << " " << s;
            }
            amrex::Print() << " (" << warp_size << " is used)\n"
                           << "  maxMemAllocSize: " << device_prop.maxMemAllocSize << "\n"
                           << "  managedMemory: " << (device_prop.managedMemory ? "Yes" : "No") << "\n"
                           << "  concurrentManagedAccess: " << (device_prop.concurrentManagedAccess ? "Yes" : "No") << "\n"
                           << "  maxParameterSize: " << device_prop.maxParameterSize << "\n"
                           << std::endl;
        }
        auto found = std::find(sgss.begin(), sgss.end(), static_cast<decltype(sgss)::value_type>(warp_size));
        if (found == sgss.end()) amrex::Abort("Incorrect subgroup size");
    }
#endif

    gpu_stream.resize(OpenMP::get_max_threads(), gpu_default_stream);

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

#ifdef AMREX_USE_DPCPP
    max_blocks_per_launch = 1000000; // xxxxx DPCPP todo
#else
    max_blocks_per_launch = numMultiProcessors() * maxThreadsPerMultiProcessor() / AMREX_GPU_MAX_THREADS;
#endif

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
    amrex::ignore_unused(idx);
#ifdef AMREX_USE_GPU
    if (idx < 0) {
        gpu_stream[OpenMP::get_thread_num()] = gpu_default_stream;

#ifdef AMREX_USE_ACC
        amrex_set_acc_stream(acc_async_sync);
#endif
    } else {
        gpu_stream[OpenMP::get_thread_num()] = gpu_streams[idx % max_gpu_streams];

#ifdef AMREX_USE_ACC
        amrex_set_acc_stream(idx % max_gpu_streams);
#endif
    }
#endif
}

#ifdef AMREX_USE_GPU
gpuStream_t
Device::resetStream () noexcept
{
    gpuStream_t r = gpu_stream[OpenMP::get_thread_num()];
    gpu_stream[OpenMP::get_thread_num()] = gpu_default_stream;
    return r;
}

gpuStream_t
Device::setStream (gpuStream_t s) noexcept
{
    gpuStream_t r = gpu_stream[OpenMP::get_thread_num()];
    gpu_stream[OpenMP::get_thread_num()] = s;
    return r;
}
#endif

void
Device::synchronize () noexcept
{
#ifdef AMREX_USE_DPCPP
    nonNullStreamSynchronize();
    try {
        gpu_default_stream.queue->wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("synchronize: ")+ex.what()+"!!!!!");
    }
#else
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipDeviceSynchronize());,
                       AMREX_CUDA_SAFE_CALL(cudaDeviceSynchronize()); )
#endif
}

void
Device::streamSynchronize () noexcept
{
#ifdef AMREX_USE_DPCPP
    auto& q = streamQueue();
    try {
        q.wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("streamSynchronize: ")+ex.what()+"!!!!!");
    }
#else
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipStreamSynchronize(gpu_stream[OpenMP::get_thread_num()]));,
                       AMREX_CUDA_SAFE_CALL(cudaStreamSynchronize(gpu_stream[OpenMP::get_thread_num()])); )
#endif
}

#ifdef AMREX_USE_DPCPP
void
Device::nonNullStreamSynchronize () noexcept
{
    for (auto const& s : gpu_streams) {
        try {
            s.queue->wait_and_throw();
        } catch (sycl::exception const& ex) {
            amrex::Abort(std::string("nonNullStreamSynchronize: ")+ex.what()+"!!!!!");
        }
    }
}
#endif

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
    amrex::ignore_unused(p,sz,device);
    // HIP does not support memory advise.
#ifdef AMREX_USE_CUDA
#ifndef AMREX_USE_HIP
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
#endif
    {
        AMREX_CUDA_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device));
    }
#elif defined(AMREX_USE_DPCPP)
    // xxxxx DPCPP todo: mem_advise
    // if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    // {
    //     auto& q = Gpu::Device::streamQueue();
    //     q.mem_advise(p, sz, PI_MEM_ADVICE_SET_PREFERRED_LOCATION);
    // }
#endif
}

void
Device::mem_advise_set_readonly (void* p, const std::size_t sz)
{
    amrex::ignore_unused(p,sz);
    // HIP does not support memory advise.
#ifdef AMREX_USE_CUDA
#ifndef AMREX_USE_HIP
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
#endif
    {
        AMREX_CUDA_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId));
    }
#elif defined(AMREX_USE_DPCPP)
    // xxxxx DPCPP todo: mem_advise
    // if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    // {
    //     auto& q = Gpu::Device::streamQueue();
    //     q.mem_advise(p, sz, PI_MEM_ADVICE_SET_READ_MOSTLY);
    // }
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
Device::n_threads_and_blocks (const Long N, dim3& numBlocks, dim3& numThreads) noexcept
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

        // Default to only an x loop in most cases.
        if (numThreadsMin.y == 1 && numThreadsMin.z == 1) {

            numBlocks.x = SM_mult_factor * num_SMs;
            numBlocks.y = 1;
            numBlocks.z = 1;

        } else {

            numBlocks.x = 1;
            numBlocks.y = SM_mult_factor;
            numBlocks.z = num_SMs;

        }

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

        // Default to only an x loop in most cases.
        if (numThreadsMin.y == 1 && numThreadsMin.z == 1) {

            numBlocks.x = SM_mult_factor * num_SMs;
            numBlocks.y = 1;
            numBlocks.z = 1;

        } else {

            numBlocks.x = 1;
            numBlocks.y = SM_mult_factor;
            numBlocks.z = num_SMs;

        }

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
    AMREX_HIP_OR_CUDA_OR_DPCPP( AMREX_HIP_SAFE_CALL(hipMemGetInfo(&f,&t));,
                                AMREX_CUDA_SAFE_CALL(cudaMemGetInfo(&f,&t));,
                                f = device_prop.totalGlobalMem; ); // xxxxx DPCPP todo
    amrex::ignore_unused(t);
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
