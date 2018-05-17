
#include <iostream>
#include <map>
#include <AMReX_Device.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

bool amrex::Device::in_device_launch_region = false;
int amrex::Device::device_id = 0;

int amrex::Device::verbose = 0;

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
cudaStream_t amrex::Device::cuda_streams[max_cuda_streams];
cudaStream_t amrex::Device::cuda_stream;

dim3 amrex::Device::numThreadsMin = dim3(1, 1, 1);

dim3 amrex::Device::numThreadsOverride = dim3(1, 1, 1);
dim3 amrex::Device::numBlocksOverride = dim3(1, 1, 1);
#endif

void
amrex::Device::initialize_cuda_c () {

    for (int i = 0; i < max_cuda_streams; ++i)
        CudaAPICheck(cudaStreamCreate(&cuda_streams[i]));

    cuda_stream = cuda_streams[0];

    ParmParse pp("device");

    int nx = 1;
    int ny = 1;
    int nz = 1;

    pp.query("numThreads.x", nx);
    pp.query("numThreads.y", ny);
    pp.query("numThreads.z", nz);

    numThreadsOverride.x = (int) nx;
    numThreadsOverride.y = (int) ny;
    numThreadsOverride.z = (int) nz;

    nx = 1;
    ny = 1;
    nz = 1;

    pp.query("numBlocks.x", nx);
    pp.query("numBlocks.y", ny);
    pp.query("numBlocks.z", nz);

    numBlocksOverride.x = (int) nx;
    numBlocksOverride.y = (int) ny;
    numBlocksOverride.z = (int) nz;

}

cudaStream_t
amrex::Device::stream_from_index(int idx) {

    if (idx < 0)
        return 0;
    else
        return cuda_streams[idx % max_cuda_streams];

}


void
amrex::Device::initialize_device() {

    ParmParse pp("device");

    pp.query("v", verbose);
    pp.query("verbose", verbose);

#ifdef AMREX_USE_CUDA

    const int n_procs = ParallelDescriptor::NProcs();
    const int my_rank = ParallelDescriptor::MyProc();
    const int ioproc  = ParallelDescriptor::IOProcessorNumber();

#ifdef AMREX_USE_NVML

    // Temporary character buffer for CUDA/NVML APIs.
    unsigned int char_length = 256;
    char c[char_length];

    NvmlAPICheck(nvmlInit());

    // Count the number of NVML visible devices.

    unsigned int nvml_device_count;
    NvmlAPICheck(nvmlDeviceGetCount(&nvml_device_count));

    if (nvml_device_count <= 0)
        return;

    // Get the PCI bus ID and UUID for each NVML visible device.

    std::vector<std::string> nvml_PCI_ids;
    std::vector<std::string> nvml_UUIDs;

    for (int i = 0; i < nvml_device_count; ++i) {

	nvmlDevice_t handle;
	NvmlAPICheck(nvmlDeviceGetHandleByIndex(i, &handle));

	NvmlAPICheck(nvmlDeviceGetUUID(handle, c, char_length));
	nvml_UUIDs.push_back(std::string(c));

        nvmlPciInfo_t pci_info;
        NvmlAPICheck(nvmlDeviceGetPciInfo(handle, &pci_info));
        nvml_PCI_ids.push_back(std::string(pci_info.busId));

    }

    // Count the number of CUDA visible devices.

    int cuda_device_count;
    CudaAPICheck(cudaGetDeviceCount(&cuda_device_count));

    if (cuda_device_count <= 0)
        return;

    // Get the PCI bus ID for each CUDA visible device.

    std::vector<std::string> cuda_PCI_ids;

    for (int cuda_device = 0; cuda_device < cuda_device_count; ++cuda_device) {

        CudaAPICheck(cudaDeviceGetPCIBusId(c, char_length, cuda_device));

        // Reset the device after accessing it; this is necessary
        // if the device is using exclusive process mode (without MPS).

        CudaAPICheck(cudaDeviceReset());

        cuda_PCI_ids.push_back(std::string(c));

    }

    // Using the PCI bus ID as a translation factor,
    // figure out the UUIDs of the CUDA visible devices.

    std::vector<std::string> cuda_UUIDs;

    for (unsigned int nvml_device = 0; nvml_device < nvml_device_count; ++nvml_device) {

        for (int cuda_device = 0; cuda_device < cuda_device_count; ++cuda_device) {

            if (cuda_PCI_ids[cuda_device] == nvml_PCI_ids[nvml_device])
                cuda_UUIDs.push_back(nvml_UUIDs[nvml_device]);

        }

    }

    // Gather the list of device UUID's from all ranks. For simplicitly we'll
    // assume that every rank sees the same number of devices and every device
    // UUID has the same number of characters; this assumption could be lifted
    // later in situations with unequal distributions of devices per node/rank.

    int strlen = cuda_UUIDs[0].size();
    int len = cuda_device_count * strlen;

    char sendbuf[len];
    char recvbuf[n_procs][len];

    int i = 0;
    for (std::string s : cuda_UUIDs) {
        for (int j = 0; j < strlen; ++j) {
            sendbuf[i * strlen + j] = s[j];
        }
        ++i;
    }

#ifdef BL_USE_MPI
    MPI_Allgather(sendbuf, len, MPI_CHAR, recvbuf, len, MPI_CHAR, MPI_COMM_WORLD);
#else
    for (int j = 0; j < len; ++j) {
        recvbuf[0][j] = sendbuf[j];
    }
#endif

    // Count up the number of ranks that share each GPU, and record their rank number.

    std::map<std::string, std::vector<int>> uuid_to_rank;

    for (std::string s : cuda_UUIDs) {
        for (int i = 0; i < n_procs; ++i) {
            for (int j = 0; j < cuda_device_count; ++j) {
                std::string temp_s(&recvbuf[i][j * strlen], strlen);
                if (s == temp_s) {
                    uuid_to_rank[s].push_back(i);
                }
            }
        }
    }

    // For each rank that shares a GPU, use round-robin assignment
    // to assign MPI ranks to GPUs. We will arbitrarily assign
    // ranks to GPUs. It would be nice to do better here and be
    // socket-aware, but this is complicated to get right. It is
    // better for this to be offloaded to the job launcher. Note that
    // the logic here will still work even if there's only one GPU
    // visible to each rank.

    device_id = -1;
    int j = 0;
    for (auto it = uuid_to_rank.begin(); it != uuid_to_rank.end(); ++it) {
        for (int i = 0; i < it->second.size(); ++i) {
            if (it->second[i] == my_rank && i % cuda_device_count == j) {

                device_id = j;

                break;
            }
        }
        j += 1;
        if (device_id != -1) break;
    }
    if (device_id == -1 || device_id >= cuda_device_count)
        amrex::Abort("Could not associate MPI rank with GPU");

    ParallelDescriptor::Barrier();

    // Total number of GPUs seen by all ranks.

    Real count = 0.0;
    for (auto it = uuid_to_rank.begin(); it != uuid_to_rank.end(); ++it) {
        count += 1.0 / it->second.size();
    }

    ParallelDescriptor::ReduceRealSum(count);

    int total_count = (int) count;

#else

    // If we don't have NVML, assign every processor to device 0.
    // The user will be on their own for ensuring that their process
    // only sees the GPU it wants in CUDA_VISIBLE_DEVICES.

    int cuda_device_count;
    CudaAPICheck(cudaGetDeviceCount(&cuda_device_count));

    if (cuda_device_count <= 0)
        return;

    device_id = 0;

    int total_count = n_procs;

#endif

    initialize_cuda(&device_id, &my_rank, &total_count, &n_procs, &ioproc, &verbose);

    initialize_cuda_c();

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
#endif

    cuda_stream = cuda_streams[idx % max_cuda_streams + 1];

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

    void* ptr = NULL;

#ifdef AMREX_USE_CUDA
    gpu_malloc(&ptr, &sz);
#endif

    return ptr;

}

void
amrex::Device::device_free(void* ptr) {

#ifdef AMREX_USE_CUDA
    gpu_free(ptr);
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
#ifndef NO_CUDA_8
    CudaAPICheck(cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device));
#endif
#endif

}

void
amrex::Device::mem_advise_set_readonly(void* p, const std::size_t sz) {

#ifdef AMReX_USE_CUDA
#ifndef NO_CUDA_8
    CudaAPICheck(cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId));
#endif
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

    numThreads.x = std::max(numThreadsMin.x, CUDA_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::max(numThreadsMin.y, CUDA_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::max(numThreadsMin.z, CUDA_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    // Allow the user to override these at runtime.

    numBlocks.x = std::max(numBlocks.x, numBlocksOverride.x);
    numBlocks.y = std::max(numBlocks.y, numBlocksOverride.y);
    numBlocks.z = std::max(numBlocks.z, numBlocksOverride.z);

    numThreads.x = std::max(numThreads.x, numThreadsOverride.x);
    numThreads.y = std::max(numThreads.y, numThreadsOverride.y);
    numThreads.z = std::max(numThreads.z, numThreadsOverride.z);

}
#endif

