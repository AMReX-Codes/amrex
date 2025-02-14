
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Machine.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#ifdef AMREX_USE_HYPRE
#  include <_hypre_utilities.h>
#endif

#include <iostream>
#include <map>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <exception>

#if defined(AMREX_USE_CUDA)
#include <cuda_profiler_api.h>
#if defined(AMREX_PROFILING) || defined (AMREX_TINY_PROFILING)
#if __has_include(<nvtx3/nvToolsExt.h>)
#  include <nvtx3/nvToolsExt.h>
#else
#  include <nvToolsExt.h>
#endif
#endif
#endif

#if defined(AMREX_USE_HIP)
#include <hip/hip_runtime.h>
#if defined(AMREX_USE_ROCTX)
#include <roctracer/roctracer_ext.h>
#if defined(AMREX_PROFILING) || defined (AMREX_TINY_PROFILING)
#include <roctracer/roctx.h>
#endif
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

#ifdef AMREX_USE_SYCL
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

#ifdef AMREX_USE_HIP
namespace {
    __host__ __device__ void amrex_check_wavefront_size () {
#ifdef __HIP_DEVICE_COMPILE__
        // * https://github.com/AMReX-Codes/amrex/issues/3792
        //   __AMDGCN_WAVEFRONT_SIZE is valid in device code only.
        //   Thus we have to check it this way.
        // * https://github.com/AMReX-Codes/amrex/issues/4270
        //   __AMDGCN_WAVEFRONT_SIZE will be deprecated.
#ifdef __AMDGCN_WAVEFRONT_SIZE
        static_assert(__AMDGCN_WAVEFRONT_SIZE == AMREX_AMDGCN_WAVEFRONT_SIZE,
                      "Please let the amrex team know if you encounter this");
#endif
#endif
    }
}
#endif

namespace amrex::Gpu {

int Device::device_id = 0;
int Device::num_devices_used = 0;
int Device::num_device_partners = 1;
int Device::verbose = 0;
#ifdef AMREX_USE_GPU
int Device::max_gpu_streams = 4;
#else
int Device::max_gpu_streams = 1;
#endif

#ifdef AMREX_USE_GPU
dim3 Device::numThreadsMin      = dim3(1, 1, 1);
dim3 Device::numThreadsOverride = dim3(0, 0, 0);
dim3 Device::numBlocksOverride  = dim3(0, 0, 0);
unsigned int Device::max_blocks_per_launch = 2560;

Vector<gpuStream_t> Device::gpu_stream_pool;
Vector<gpuStream_t> Device::gpu_stream;
gpuDeviceProp_t     Device::device_prop;
int                 Device::memory_pools_supported = 0;

constexpr int Device::warp_size;

#ifdef AMREX_USE_SYCL
std::unique_ptr<sycl::context> Device::sycl_context;
std::unique_ptr<sycl::device>  Device::sycl_device;
#endif

namespace {

#if defined(__CUDACC__)
    AMREX_GPU_GLOBAL void emptyKernel() {}
#endif

    void InitializeGraph(int graph_size)
    {
        amrex::ignore_unused(graph_size);

#if defined(__CUDACC__) && defined(AMREX_USE_CUDA)

        BL_PROFILE("InitGraph");

        int streams = Gpu::Device::numGpuStreams();
        cudaGraphExec_t graphExec{};
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
#ifdef AMREX_USE_GPU

#if defined(AMREX_USE_CUDA) && (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
    // Wrap cuda init to identify it appropriately in nvvp.
    // Note: first substantial cuda call may cause a lengthy
    // cuda API and cuda driver API initialization that will
    // be captured by the profiler. It a necessary, system
    // dependent step that is unavoidable.
    nvtxRangePush("initialize_device");
#endif

    ParmParse ppamrex("amrex");
    ppamrex.queryAdd("max_gpu_streams", max_gpu_streams);
    max_gpu_streams = std::min(max_gpu_streams, AMREX_GPU_MAX_STREAMS);
    max_gpu_streams = std::max(max_gpu_streams, 1);

    ParmParse pp("device");
    if (! pp.query("verbose", "v", verbose)) {
        pp.add("verbose", verbose);
    }

    if (amrex::Verbose()) {
        AMREX_HIP_OR_CUDA_OR_SYCL
            ( amrex::Print() << "Initializing HIP...\n";,
              amrex::Print() << "Initializing CUDA...\n";,
              amrex::Print() << "Initializing SYCL...\n"; )
    }

    // Count the number of GPU devices.
    int gpu_device_count = 0;
#ifdef AMREX_USE_SYCL
    {
        sycl::platform platform(sycl::gpu_selector_v);
        auto const& gpu_devices = platform.get_devices();
        gpu_device_count = gpu_devices.size();
        if (gpu_device_count <= 0) {
            amrex::Abort("No GPU device found");
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

    int n_local_procs = 1;
    amrex::ignore_unused(n_local_procs);

    if (ParallelDescriptor::NProcs() == 1) {
        device_id = 0;
    }
    else if (gpu_device_count == 1) {
        device_id = 0;
    }
    else {
        if (ParallelDescriptor::NProcsPerNode() == gpu_device_count) {
            device_id = ParallelDescriptor::MyRankInNode();
        } else if (ParallelDescriptor::NProcsPerProcessor() == gpu_device_count) {
            device_id = ParallelDescriptor::MyRankInProcessor();
        } else {
            device_id = ParallelDescriptor::MyProc() % gpu_device_count;
        }
    }

    if (gpu_device_count > 1){
        if (Machine::name() == "nersc.perlmutter") {
            // The CPU/GPU mapping on perlmutter has the reverse order.
            device_id = gpu_device_count - device_id - 1;
            if (amrex::Verbose()) {
                amrex::Print() << "Multiple GPUs are visible to each MPI rank. Fixing GPU assignment for Perlmuuter according to heuristics.\n";
            }
        } else if (Machine::name() == "olcf.frontier") {
            // The CPU/GPU mapping on fronter is documented at
            // https://docs.olcf.ornl.gov/systems/frontier_user_guide.html
            if (gpu_device_count == 8) {
                constexpr std::array<int,8> gpu_order = {4,5,2,3,6,7,0,1};
                device_id = gpu_order[device_id];
                if (amrex::Verbose()) {
                    amrex::Print() << "Multiple GPUs are visible to each MPI rank. Fixing GPU assignment for Frontier according to heuristics.\n";
                }
            }
        } else {
            if (amrex::Verbose() && ParallelDescriptor::IOProcessor()) {
                amrex::Warning("Multiple GPUs are visible to each MPI rank. This is usually not an issue. But this may lead to incorrect or suboptimal rank-to-GPU mapping.");
            }
        }
    }

    AMREX_HIP_OR_CUDA(AMREX_HIP_SAFE_CALL (hipSetDevice(device_id));,
                      AMREX_CUDA_SAFE_CALL(cudaSetDevice(device_id)); );

#ifdef AMREX_USE_ACC
    amrex_initialize_acc(device_id);
#endif

    initialize_gpu();

    num_devices_used = ParallelDescriptor::NProcs();

#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {

#if defined(HIP_VERSION_MAJOR) && defined(HIP_VERSION_MINOR) && ((HIP_VERSION_MAJOR < 5) || ((HIP_VERSION_MAJOR == 5) && (HIP_VERSION_MINOR < 2)))

        // hip < 5.2: uuid not supported
        num_device_partners = 1;

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

        constexpr int len = 16;
        static_assert(std::is_same<decltype(AMREX_HIP_OR_CUDA(hipUUID,cudaUUID_t)::bytes),
                                   char[len]>());
        std::vector<char> buf(ParallelDescriptor::NProcs()*len);
        char* pbuf = buf.data();
#ifdef AMREX_USE_CUDA
        auto const& uuid = device_prop.uuid;
#else
        hipUUID uuid;
        AMREX_HIP_SAFE_CALL(hipDeviceGetUuid(&uuid, device_id));
#endif
        char const* sbuf = uuid.bytes;
        MPI_Allgather(sbuf, len, MPI_CHAR, pbuf, len, MPI_CHAR,
                      ParallelDescriptor::Communicator());
        std::map<std::string,int> uuid_counts;
        std::string my_uuid;
        for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
            std::string iuuid(pbuf+i*len, len);
            if (i == ParallelDescriptor::MyProc()) {
                my_uuid = iuuid;
            }
            ++uuid_counts[iuuid];
        }
        num_devices_used = uuid_counts.size();
        num_device_partners = uuid_counts[my_uuid];

#elif defined(AMREX_USE_SYCL)

        auto const& d = *sycl_device;
        if (d.has(sycl::aspect::ext_intel_device_info_uuid)) {
            auto uuid = d.get_info<sycl::ext::intel::info::device::uuid>();
            using id_t = decltype(uuid); // std::array<unsigned char,16>
            using char_t = id_t::value_type; // unsigned char
            int len = std::tuple_size<id_t>::value;
            std::vector<char_t> buf(ParallelDescriptor::NProcs()*len);
            char_t* pbuf = buf.data();
            MPI_Allgather(uuid.data(), len,
                          ParallelDescriptor::Mpi_typemap<char_t>::type(),
                          pbuf, len,
                          ParallelDescriptor::Mpi_typemap<char_t>::type(),
                          ParallelDescriptor::Communicator());
            using str_t = std::basic_string<char_t>;
            std::map<str_t,int> uuid_counts;
            str_t my_uuid;
            for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
                str_t iuuid(pbuf+i*len, len);
                if (i == ParallelDescriptor::MyProc()) {
                    my_uuid = iuuid;
                }
                ++uuid_counts[iuuid];
            }
            num_devices_used = uuid_counts.size();
            num_device_partners = uuid_counts[my_uuid];

#if 0
            for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
                if (i == ParallelDescriptor::MyProc()) {
                    std::cout << "Proc. " << i << ": |";
                    for (auto x : uuid) {
                        std::cout << std::hex << static_cast<unsigned int>(x) << "|";
                    }
                    std::cout << '\n';;
                }
                ParallelDescriptor::Barrier();
            }
#endif
        }

#endif

        AMREX_ALWAYS_ASSERT(num_device_partners > 0);
    }
#endif /* AMREX_USE_MPI */

    if (amrex::Verbose()) {
#if defined(AMREX_USE_CUDA)
        amrex::Print() << "CUDA"
#elif defined(AMREX_USE_HIP)
        amrex::Print() << "HIP"
#elif defined(AMREX_USE_SYCL)
        amrex::Print() << "SYCL"
#endif
                       << " initialized with " << num_devices_used
                       << ((num_devices_used == 1) ? " device.\n"
                                                   : " devices.\n");
        if (num_devices_used < ParallelDescriptor::NProcs() && ParallelDescriptor::IOProcessor()) {
            amrex::Warning("There are more MPI processes than the number of unique GPU devices. This is not necessarily a problem.\n"
                           "For example this could happen when a device such as MI300A is partitioned into multiple subdevices.");
        }
    }

#if defined(AMREX_USE_CUDA) && (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
    nvtxRangePop();
#endif

#if defined(AMREX_USE_HIP)
    if (num_devices_used < 0) {
        // This test is always false, but it makes the compiler no longer
        // complain about unused function, amrex_check_wavefront_size.
        amrex::single_task(amrex_check_wavefront_size);
    }
#endif

    Device::profilerStart();

#endif /* AMREX_USE_GPU */
}

void
Device::Finalize ()
{
#ifdef AMREX_USE_GPU
    Device::profilerStop();

#ifdef AMREX_USE_SYCL
    for (auto& s : gpu_stream_pool) {
        delete s.queue;
        s.queue = nullptr;
    }
    sycl_context.reset();
    sycl_device.reset();
#else
    for (int i = 0; i < max_gpu_streams; ++i)
    {
        AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL( hipStreamDestroy(gpu_stream_pool[i]));,
                          AMREX_CUDA_SAFE_CALL(cudaStreamDestroy(gpu_stream_pool[i])); );
    }
#endif

    gpu_stream.clear();

#ifdef AMREX_USE_ACC
    amrex_finalize_acc();
#endif

#endif
}

void
Device::initialize_gpu ()
{
#ifdef AMREX_USE_GPU

    gpu_stream_pool.resize(max_gpu_streams);

#ifdef AMREX_USE_HIP

    AMREX_HIP_SAFE_CALL(hipGetDeviceProperties(&device_prop, device_id));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warp_size == device_prop.warpSize, "Incorrect warp size");

    // check compute capability

    // AMD devices do not support shared cache banking.

    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_HIP_SAFE_CALL(hipStreamCreate(&gpu_stream_pool[i]));
    }

#ifdef AMREX_GPU_STREAM_ALLOC_SUPPORT
    AMREX_HIP_SAFE_CALL(hipDeviceGetAttribute(&memory_pools_supported, hipDeviceAttributeMemoryPoolsSupported, device_id));
#endif

#elif defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(device_prop.major >= 4 || (device_prop.major == 3 && device_prop.minor >= 5),
                                     "Compute capability must be >= 3.5");

#ifdef AMREX_GPU_STREAM_ALLOC_SUPPORT
    cudaDeviceGetAttribute(&memory_pools_supported, cudaDevAttrMemoryPoolsSupported, device_id);
#endif

#if (__CUDACC_VER_MAJOR__ < 12) || ((__CUDACC_VER_MAJOR__ == 12) && (__CUDACC_VER_MINOR__ < 4))
    if (sizeof(Real) == 8) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    } else if (sizeof(Real) == 4) {
        AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
    }
#endif

    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_CUDA_SAFE_CALL(cudaStreamCreate(&gpu_stream_pool[i]));
#ifdef AMREX_USE_ACC
        acc_set_cuda_stream(i, gpu_stream_pool[i]);
#endif
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warp_size == device_prop.warpSize, "Incorrect warp size");

#elif defined(AMREX_USE_SYCL)
    { // create device, context and queues
        sycl::platform platform(sycl::gpu_selector_v);
        auto const& gpu_devices = platform.get_devices();
        sycl_device = std::make_unique<sycl::device>(gpu_devices[device_id]);
        sycl_context = std::make_unique<sycl::context>(*sycl_device, amrex_sycl_error_handler);
        for (int i = 0; i < max_gpu_streams; ++i) {
            gpu_stream_pool[i].queue = new sycl::queue(*sycl_context, *sycl_device,
                                         sycl::property_list{sycl::property::queue::in_order{}});
        }
    }

    { // device property
        auto const& d = *sycl_device;
        device_prop.name = d.get_info<sycl::info::device::name>();
        device_prop.vendor = d.get_info<sycl::info::device::vendor>();
        device_prop.totalGlobalMem = d.get_info<sycl::info::device::global_mem_size>();
        device_prop.sharedMemPerBlock = d.get_info<sycl::info::device::local_mem_size>();
        device_prop.multiProcessorCount = d.get_info<sycl::info::device::max_compute_units>();
        device_prop.maxThreadsPerMultiProcessor = -1; // xxxxx SYCL todo: d.get_info<sycl::info::device::max_work_items_per_compute_unit>(); // unknown
        device_prop.maxThreadsPerBlock = d.get_info<sycl::info::device::max_work_group_size>();
        auto mtd = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
        device_prop.maxThreadsDim[0] = mtd[0];
        device_prop.maxThreadsDim[1] = mtd[1];
        device_prop.maxThreadsDim[2] = mtd[2];
        device_prop.maxGridSize[0] = -1; // xxxxx SYCL todo: unknown
        device_prop.maxGridSize[0] = -1; // unknown
        device_prop.maxGridSize[0] = -1; // unknown
        device_prop.warpSize = warp_size;
        auto sgss = d.get_info<sycl::info::device::sub_group_sizes>();
        device_prop.maxMemAllocSize = d.get_info<sycl::info::device::max_mem_alloc_size>();
        device_prop.managedMemory = d.has(sycl::aspect::usm_host_allocations);
        device_prop.concurrentManagedAccess = d.has(sycl::aspect::usm_shared_allocations);
        device_prop.maxParameterSize = d.get_info<sycl::info::device::max_parameter_size>();
        if (verbose)
        {
            amrex::Print() << "Device Properties:\n"
                           << "  name: " << device_prop.name << "\n"
                           << "  vendor: " << device_prop.vendor << "\n"
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
                           << '\n';
#if defined(__INTEL_LLVM_COMPILER)
            if (d.has(sycl::aspect::ext_intel_gpu_eu_simd_width)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_eu_simd_width>();
                amrex::Print() << "  Intel GPU Execution Unit SIMD Width: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_gpu_eu_count)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_eu_count>();
                amrex::Print() << "  Intel GPU Execution Unit Count: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_gpu_slices)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_slices>();
                amrex::Print() << "  Intel GPU Number of Slices: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_gpu_subslices_per_slice)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_subslices_per_slice>();
                amrex::Print() << "  Intel GPU Number of Subslices per Slice: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_gpu_eu_count_per_subslice)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_eu_count_per_subslice>();
                amrex::Print() << "  Intel GPU Number of EUs per Subslice: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_gpu_hw_threads_per_eu)) {
                auto r = d.get_info<sycl::ext::intel::info::device::gpu_hw_threads_per_eu>();
                amrex::Print() << "  Intel GPU Number of hardware threads per EU: " << r << "\n";
            }
            if (d.has(sycl::aspect::ext_intel_max_mem_bandwidth)) {
                auto r = d.get_info<sycl::ext::intel::info::device::max_mem_bandwidth>();
                amrex::Print() << "  Intel GPU Maximum Memory Bandwidth (B/s): " << r << "\n";
            }
#endif
        }
        auto found = std::find(sgss.begin(), sgss.end(), static_cast<decltype(sgss)::value_type>(warp_size));
        if (found == sgss.end()) { amrex::Abort("Incorrect subgroup size"); }
    }
#endif

    gpu_stream.resize(OpenMP::get_max_threads(), gpu_stream_pool[0]);

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

#ifdef AMREX_USE_SYCL
    // max_blocks_per_launch = 100000; // xxxxx SYCL todo
#else
    max_blocks_per_launch = 4 * numMultiProcessors() * maxThreadsPerMultiProcessor() / AMREX_GPU_MAX_THREADS;
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

int Device::numDevicePartners () noexcept
{
    return num_device_partners;
}

#ifdef AMREX_USE_GPU
int
Device::streamIndex (gpuStream_t s) noexcept
{
    auto it = std::find(std::begin(gpu_stream_pool), std::end(gpu_stream_pool), s);
    return static_cast<int>(std::distance(std::begin(gpu_stream_pool), it));
}
#endif

void
Device::setStreamIndex (int idx) noexcept
{
    amrex::ignore_unused(idx);
#ifdef AMREX_USE_GPU
    gpu_stream[OpenMP::get_thread_num()] = gpu_stream_pool[idx % max_gpu_streams];
#ifdef AMREX_USE_ACC
    amrex_set_acc_stream(idx % max_gpu_streams);
#endif
#endif
}

#ifdef AMREX_USE_GPU
gpuStream_t
Device::resetStream () noexcept
{
    gpuStream_t r = gpu_stream[OpenMP::get_thread_num()];
    gpu_stream[OpenMP::get_thread_num()] = gpu_stream_pool[0];
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
#ifdef AMREX_USE_SYCL
    for (auto const& s : gpu_stream_pool) {
        try {
            s.queue->wait_and_throw();
        } catch (sycl::exception const& ex) {
            amrex::Abort(std::string("synchronize: ")+ex.what()+"!!!!!");
        }
    }
#else
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipDeviceSynchronize());,
                       AMREX_CUDA_SAFE_CALL(cudaDeviceSynchronize()); )
#endif
}

void
Device::streamSynchronize () noexcept
{
#ifdef AMREX_USE_SYCL
    auto& q = streamQueue();
    try {
        q.wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("streamSynchronize: ")+ex.what()+"!!!!!");
    }
#else
    AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipStreamSynchronize(gpuStream()));,
                       AMREX_CUDA_SAFE_CALL(cudaStreamSynchronize(gpuStream())); )
#endif
}

void
Device::streamSynchronizeAll () noexcept
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    Device::synchronize();
#else
    for (auto const& s : gpu_stream_pool) {
        AMREX_HIP_OR_CUDA( AMREX_HIP_SAFE_CALL(hipStreamSynchronize(s));,
                           AMREX_CUDA_SAFE_CALL(cudaStreamSynchronize(s)); )
    }
#endif
#endif
}

#if defined(__CUDACC__) && defined(AMREX_USE_CUDA)

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
    cudaGraphExec_t graphExec{};

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
//  Implements cudaGraphInstantiate error logging feature.
//  Upon error, delays abort until message is output.
    constexpr int log_size = 1028;
    char graph_log[log_size];
    graph_log[0]='\0';

    cudaGraphInstantiate(&graphExec, graph, NULL, &(graph_log[0]), log_size);

    if (graph_log[0] != '\0')
    {
        amrex::Print() << graph_log << '\n';
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
Device::mem_advise_set_preferred (void* p, std::size_t sz, int device)
{
    amrex::ignore_unused(p,sz,device);
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    {
        AMREX_HIP_OR_CUDA
            (AMREX_HIP_SAFE_CALL(
                 hipMemAdvise(p, sz, hipMemAdviseSetPreferredLocation, device)),
             AMREX_CUDA_SAFE_CALL(
                 cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device)));
    }
#elif defined(AMREX_USE_SYCL)
    // xxxxx SYCL todo: mem_advise
    // if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    // {
    //     auto& q = Gpu::Device::streamQueue();
    //     q.mem_advise(p, sz, PI_MEM_ADVICE_SET_PREFERRED_LOCATION);
    // }
#endif
}

void
Device::mem_advise_set_readonly (void* p, std::size_t sz)
{
    amrex::ignore_unused(p,sz);
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    {
        AMREX_HIP_OR_CUDA
            (AMREX_HIP_SAFE_CALL(
                 hipMemAdvise(p, sz, hipMemAdviseSetReadMostly, hipCpuDeviceId)),
             AMREX_CUDA_SAFE_CALL(
                 cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId)));
    }
#elif defined(AMREX_USE_SYCL)
    // xxxxx SYCL todo: mem_advise
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
    numBlocks = std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)); // in case N = 0
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

    AMREX_ASSERT(numThreads.x <= static_cast<unsigned>(device_prop.maxThreadsDim[0]));
    AMREX_ASSERT(numThreads.y <= static_cast<unsigned>(device_prop.maxThreadsDim[1]));
    AMREX_ASSERT(numThreads.z <= static_cast<unsigned>(device_prop.maxThreadsDim[2]));
    AMREX_ASSERT(numThreads.x*numThreads.y*numThreads.z <= static_cast<unsigned>(device_prop.maxThreadsPerBlock));
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

#endif

std::size_t
Device::freeMemAvailable ()
{
#ifdef AMREX_USE_GPU
    std::size_t f;
#ifndef AMREX_USE_SYCL
    std::size_t t;
#endif
    AMREX_HIP_OR_CUDA_OR_SYCL( AMREX_HIP_SAFE_CALL(hipMemGetInfo(&f,&t));,
                             AMREX_CUDA_SAFE_CALL(cudaMemGetInfo(&f,&t));,
                               f = device_prop.totalGlobalMem; );
#if defined (AMREX_USE_SYCL) && defined(__INTEL_LLVM_COMPILER)
    if (sycl_device->has(sycl::aspect::ext_intel_free_memory)) {
        f = sycl_device->get_info<sycl::ext::intel::info::device::free_memory>();
    }
#endif
    return f;
#else
    return 0;
#endif
}

void
Device::profilerStart ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStart());
#elif (defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX))
    roctracer_start();
#endif

}

void
Device::profilerStop ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStop());
#elif (defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX))
    roctracer_stop();
#endif
}

#ifdef AMREX_USE_HYPRE
void hypreSynchronize ()
{
#ifdef AMREX_USE_GPU
    hypre_SyncCudaDevice(hypre_handle()); // works for non-cuda device too
#endif
}
#endif

}
