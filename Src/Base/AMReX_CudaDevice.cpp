
#include <iostream>
#include <map>
#include <algorithm>
#include <AMReX_CudaDevice.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

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
namespace Cuda {

int Device::device_id = 0;
int Device::verbose = 0;

#if defined(AMREX_USE_CUDA)
constexpr int Device::max_cuda_streams;

std::array<cudaStream_t,Device::max_cuda_streams> Device::cuda_streams;
cudaStream_t Device::cuda_stream;

dim3 Device::numThreadsMin      = dim3(1, 1, 1);
dim3 Device::numThreadsOverride = dim3(0, 0, 0);
dim3 Device::numBlocksOverride  = dim3(0, 0, 0);

cudaDeviceProp Device::device_prop;
#endif

void
Device::Initialize ()
{
#ifdef AMREX_USE_CUDA

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
    AMREX_GPU_SAFE_CALL(cudaGetDeviceCount(&cuda_device_count));

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

    AMREX_GPU_SAFE_CALL(cudaSetDevice(device_id));
    AMREX_GPU_SAFE_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));

#ifdef AMREX_USE_ACC
    amrex_initialize_acc(device_id);
#endif

    initialize_cuda();

#if defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING)
    nvtxRangeEnd(nvtx_init);
#endif

    if (amrex::Verbose()) {
#ifdef AMREX_USE_MPI
        amrex::Print() << "CUDA initialized with 1 GPU per MPI rank\n";
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
#ifdef AMREX_USE_CUDA

    cudaProfilerStop();

    for (int i = 0; i < max_cuda_streams; ++i)
    {
        AMREX_GPU_SAFE_CALL(cudaStreamDestroy(cuda_streams[i]));
    }

#ifdef AMREX_USE_ACC
    amrex_finalize_acc();
#endif

    AMREX_GPU_SAFE_CALL(cudaDeviceReset());
#endif
}

void
Device::initialize_cuda ()
{
#if defined(AMREX_USE_CUDA)
    AMREX_GPU_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));

    if (device_prop.warpSize != 32) {
        amrex::Warning("Warp size != 32");
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(device_prop.major >= 6, "Compute capability must be >= 6");

    // Prefer L1 cache to shared memory (this has no effect on GPUs with a fixed L1 cache size).
    AMREX_GPU_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    if (sizeof(Real) == 8) {
        AMREX_GPU_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    } else if (sizeof(Real) == 4) {
        AMREX_GPU_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
    }

    for (int i = 0; i < max_cuda_streams; ++i) {
        AMREX_GPU_SAFE_CALL(cudaStreamCreate(&cuda_streams[i]));
    }

    cuda_stream = cuda_streams[0];

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
#endif
}

int
Device::deviceId () noexcept
{
    return device_id;
}

void
Device::setStreamIndex (const int idx) noexcept
{
#ifdef AMREX_USE_CUDA
    if (idx < 0) {
        cuda_stream = 0;
    } else {
        cuda_stream = cuda_streams[idx % max_cuda_streams];
    }
#endif
}

void
Device::synchronize ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaDeviceSynchronize());
#endif
}

void
Device::streamSynchronize ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaStreamSynchronize(cuda_stream));
#endif
}

void
Device::htod_memcpy (void* p_d, const void* p_h, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaMemcpy(p_d, p_h, sz, cudaMemcpyHostToDevice));
#endif

}

void
Device::dtoh_memcpy (void* p_h, const void* p_d, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaMemcpy(p_h, p_d, sz, cudaMemcpyDeviceToHost));
#endif

}

void
Device::htod_memcpy_async (void* p_d, const void* p_h, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_stream));
#endif

}

void
Device::dtoh_memcpy_async (void* p_h, const void* p_d, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_stream));
#endif

}

void
Device::mem_advise_set_preferred (void* p, const std::size_t sz, const int device) {

#ifdef AMREX_USE_CUDA
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
        AMREX_GPU_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device));
#endif

}

void
Device::mem_advise_set_readonly (void* p, const std::size_t sz) {
#ifdef AMREX_USE_CUDA
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
        AMREX_GPU_SAFE_CALL(cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId));
#endif
}

#if defined(AMREX_USE_CUDA)

void
Device::setNumThreadsMin (int nx, int ny, int nz) noexcept {
    numThreadsMin.x = nx;
    numThreadsMin.y = ny;
    numThreadsMin.z = nz;
}

void
Device::n_threads_and_blocks (const long N, dim3& numBlocks, dim3& numThreads) noexcept
{
    numThreads = AMREX_CUDA_MAX_THREADS;
    numBlocks = std::max((N + AMREX_CUDA_MAX_THREADS - 1) / AMREX_CUDA_MAX_THREADS, 1L); // in case N = 0
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

    numThreads.x = std::min(tile_size[0], AMREX_CUDA_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = 1;
    numBlocks.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(tile_size[0]), AMREX_CUDA_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(tile_size[1]), AMREX_CUDA_MAX_THREADS / numThreads.x   );
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = (tile_size[1] + numThreads.y - 1) / numThreads.y;
    numBlocks.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_CUDA_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreads.y   ));

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

    numThreads.x = std::min(device_prop.maxThreadsDim[0], AMREX_CUDA_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_CUDA_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_CUDA_MAX_THREADS / numThreads.x);
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_CUDA_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreads.y   ));

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

    numThreads.x = std::min(device_prop.maxThreadsDim[0], AMREX_CUDA_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_CUDA_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_CUDA_MAX_THREADS / numThreads.x);
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_CUDA_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_CUDA_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = std::max(numThreadsMin.z, numThreads.z);

#endif

    // Limit the number of threads per block to be no larger in each dimension
    // than the sizes of the box in the corresponding dimension.

    numThreads.x = std::min(numThreads.x, static_cast<unsigned>(bx.length(0)));
    numThreads.y = std::min(numThreads.y, static_cast<unsigned>(bx.length(1)));
    numThreads.z = std::min(numThreads.z, static_cast<unsigned>(bx.length(2)));

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
#ifdef AMREX_USE_CUDA
    std::size_t f, t;
    AMREX_GPU_SAFE_CALL(cudaMemGetInfo(&f,&t));
    return f;
#else
    return 0;
#endif
}

}}

