
#include <iostream>
#include <map>
#include <algorithm>
#include <AMReX_Device.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

int amrex::Device::device_id = 0;
int amrex::Device::verbose = 0;

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
const int amrex::Device::max_cuda_streams;

cudaStream_t amrex::Device::cuda_streams[max_cuda_streams];
cudaStream_t amrex::Device::cuda_stream;

dim3 amrex::Device::numThreadsMin      = dim3(1, 1, 1);
dim3 amrex::Device::numThreadsOverride = dim3(0, 0, 0);
dim3 amrex::Device::numBlocksOverride  = dim3(0, 0, 0);

cudaDeviceProp amrex::Device::device_prop;
#endif

#if defined(AMREX_USE_CUDA)
void
amrex::Device::initialize_cuda_c () {

    for (int i = 0; i < max_cuda_streams; ++i)
        CudaAPICheck(cudaStreamCreate(&cuda_streams[i]));

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

    cudaGetDeviceProperties(&device_prop, device_id);
}

cudaStream_t
amrex::Device::stream_from_index(int idx) {

    if (idx < 0)
        return 0;
    else
        return cuda_streams[idx % max_cuda_streams];

}
#endif

void
amrex::Device::initialize_device() {

    ParmParse pp("device");

    pp.query("v", verbose);
    pp.query("verbose", verbose);

#ifdef AMREX_USE_CUDA

    // Count the number of CUDA visible devices.

    int cuda_device_count;
    CudaAPICheck(cudaGetDeviceCount(&cuda_device_count));

    if (cuda_device_count <= 0)
        return;

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

    CudaAPICheck(cudaSetDevice(device_id));

    initialize_cuda(&device_id);

    initialize_cuda_c();

    if (amrex::Verbose())
        amrex::Print() << "CUDA initialized with 1 GPU per MPI rank\n";

#endif

}

void
amrex::Device::finalize_device() {

#ifdef AMREX_USE_CUDA
    finalize_cuda();

    set_is_program_running(0);
#endif

}

int
amrex::Device::deviceId() {

    return device_id;

}

void
amrex::Device::set_stream_index(const int idx) {

#ifdef AMREX_USE_CUDA
    set_stream_idx(idx);
    cuda_stream = cuda_streams[idx % max_cuda_streams + 1];
#endif
}

int
amrex::Device::get_stream_index() {

    int index = -1;

#ifdef AMREX_USE_CUDA
    get_stream_idx(&index);
#endif

    return index;

}

void
amrex::Device::prepare_for_launch(const int* lo, const int* hi) {

    // Sets the number of threads and blocks in Fortran.
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    int txmin = numThreadsMin.x;
    int tymin = numThreadsMin.y;
    int tzmin = numThreadsMin.z;

    set_threads_and_blocks(lo, hi, &txmin, &tymin, &tzmin);
#endif

}

void*
amrex::Device::get_host_pointer(const void* ptr) {

    void* r = const_cast<void*>(ptr);
#ifdef AMREX_USE_CUDA
    gpu_host_device_ptr(&r, ptr);
#endif
    return r;

}

void
amrex::Device::check_for_errors() {

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    CudaErrorCheck();
#endif

}

void
amrex::Device::synchronize() {

#ifdef AMREX_USE_CUDA
    gpu_synchronize();
#endif

}

void
amrex::Device::stream_synchronize(const int idx) {

#ifdef AMREX_USE_CUDA
    gpu_stream_synchronize(idx);
#endif

}

void*
amrex::Device::device_malloc(const std::size_t sz) {

    void* ptr = nullptr;
#ifdef AMREX_USE_CUDA
    gpu_malloc(&ptr, &sz);
#else
    ptr = amrex_malloc(sz);
#endif

    return ptr;
}

bool
amrex::Device::checkManaged(const void* ptr) {

#ifdef AMREX_USE_CUDA
     cudaPointerAttributes ptr_attr;
     cudaPointerGetAttributes(&ptr_attr, ptr);
     cudaError_t err = cudaGetLastError(); 
     if (err == cudaErrorInvalidValue)
     {
        std::cout << " cudaErrorInvalidValue ";
        return false;
     }
     return ptr_attr.isManaged;
#else
     return false;
#endif

}

bool
amrex::Device::checkDevicePtr (const void* ptr) {

#ifdef AMREX_USE_CUDA
     cudaPointerAttributes ptr_attr;
     cudaPointerGetAttributes(&ptr_attr, ptr);
     cudaError_t err = cudaGetLastError(); 
     if (err == cudaErrorInvalidValue)
     {
        std::cout << " cudaErrorInvalidValue ";
        return false;
     }
     else if (ptr_attr.memoryType == cudaMemoryTypeHost)
     {
        std::cout << " cudaMemoryTypeHost ";
        return false;
     }

     return true;
#else
     return false;
#endif

}

void
amrex::Device::device_free(void* ptr) {

#ifdef AMREX_USE_CUDA
    gpu_free(ptr);
#else
    amrex_free(ptr);
#endif

}

void
amrex::Device::device_htod_memcpy(void* p_d, const void* p_h, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    CudaAPICheck(cudaMemcpy(p_d, p_h, sz, cudaMemcpyHostToDevice));
#endif

}

void
amrex::Device::device_dtoh_memcpy(void* p_h, const void* p_d, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    CudaAPICheck(cudaMemcpy(p_h, p_d, sz, cudaMemcpyDeviceToHost));
#endif

}

void
amrex::Device::device_htod_memcpy_async(void* p_d, const void* p_h, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    CudaAPICheck(cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_stream));
#endif

}

void
amrex::Device::device_dtoh_memcpy_async(void* p_h, const void* p_d, const std::size_t sz) {

#ifdef AMREX_USE_CUDA
    CudaAPICheck(cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_stream));
#endif

}

void
amrex::Device::mem_advise_set_preferred(void* p, const std::size_t sz, const int device) {

#ifdef AMREX_USE_CUDA
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
        CudaAPICheck(cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device));
#endif

}

void
amrex::Device::mem_advise_set_readonly(void* p, const std::size_t sz) {
#ifdef AMREX_USE_CUDA
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
        CudaAPICheck(cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId));
#endif
}

void
amrex::Device::start_profiler() {

#ifdef AMREX_USE_CUDA
    gpu_start_profiler();
#endif

}

void
amrex::Device::stop_profiler() {

#ifdef AMREX_USE_CUDA
    gpu_stop_profiler();
#endif

}

#if (defined(AMREX_USE_CUDA) && defined(__CUDACC__))
void
amrex::Device::c_comps_threads_and_blocks(const int* lo, const int* hi, const int comps, dim3& numBlocks, dim3& numThreads) {

    c_threads_and_blocks(lo, hi, numBlocks, numThreads);
    numBlocks.x *= static_cast<unsigned>(comps);
}

void
amrex::Device::c_threads_and_blocks(const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) {

    int bx, by, bz, tx, ty, tz;

    int txmin = numThreadsMin.x;
    int tymin = numThreadsMin.y;
    int tzmin = numThreadsMin.z;

    get_threads_and_blocks(lo, hi, &bx, &by, &bz, &tx, &ty, &tz, &txmin, &tymin, &tzmin);

    numBlocks.x = bx;
    numBlocks.y = by;
    numBlocks.z = bz;

    numThreads.x = tx;
    numThreads.y = ty;
    numThreads.z = tz;

}

void
amrex::Device::grid_stride_threads_and_blocks(dim3& numBlocks, dim3& numThreads) {

    int num_SMs;

    int SM_mult_factor = 32;

    get_num_SMs(&num_SMs);

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

    numThreads.x = std::max((int) numThreadsMin.x, 16);
    numThreads.y = std::max((int) numThreadsMin.y, 16);
    numThreads.z = std::max((int) numThreadsMin.z, 1);

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
amrex::Device::particle_threads_and_blocks(const int np, int& numThreads, int& numBlocks) {
    numThreads = 256;
    numBlocks = (np + 256 - 1) / 256;
}


void
amrex::Device::n_threads_and_blocks (const int N, dim3& numBlocks, dim3& numThreads)
{
    const int maxBlockSize = 256;
    numThreads = maxBlockSize;
    numBlocks = std::min((N + maxBlockSize - 1) / maxBlockSize, 1); // in case N = 0
}
#endif

extern "C" {
    void* amrex_gpu_malloc (std::size_t size)
    {
#ifdef AMREX_USE_CUDA
        void *ptr = nullptr;
        cudaMalloc((void**) ptr, size); 
        return ptr;
#else
        return amrex_malloc(size);
#endif
    }


    void amrex_gpu_free (void* p)
    {
#ifdef AMREX_USE_CUDA
        cudaFree(p); 
#else
        amrex_free(p);
#endif
    }
}
