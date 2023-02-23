
#include <AMReX_Arena.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>
#include <AMReX_PArena.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>

#ifdef _WIN32
///#include <memoryapi.h>
//#define AMREX_MLOCK(x,y) VirtualLock(x,y)
//#define AMREX_MUNLOCK(x,y) VirtualUnlock(x,y)
//#define AMREX_MLOCK(x,y) ((void)0)
#define AMREX_MUNLOCK(x,y) ((void)0)
#else
#include <sys/mman.h>
//#define AMREX_MLOCK(x,y) mlock(x,y)
#define AMREX_MUNLOCK(x,y) munlock(x,y)
#endif

namespace amrex {

namespace {
    bool initialized = false;

    Arena* the_arena = nullptr;
    Arena* the_async_arena = nullptr;
    Arena* the_device_arena = nullptr;
    Arena* the_managed_arena = nullptr;
    Arena* the_pinned_arena = nullptr;
    Arena* the_cpu_arena = nullptr;

    Long the_arena_init_size = 0L;
    Long the_device_arena_init_size = 1024*1024*8;
    Long the_managed_arena_init_size = 1024*1024*8;
    Long the_pinned_arena_init_size = 1024*1024*8;
    Long the_arena_release_threshold = std::numeric_limits<Long>::max();
    Long the_device_arena_release_threshold = std::numeric_limits<Long>::max();
    Long the_managed_arena_release_threshold = std::numeric_limits<Long>::max();
    Long the_pinned_arena_release_threshold = std::numeric_limits<Long>::max();
    Long the_async_arena_release_threshold = std::numeric_limits<Long>::max();
#ifdef AMREX_USE_HIP
    bool the_arena_is_managed = false; // xxxxx HIP FIX HERE
#else
    bool the_arena_is_managed = true;
#endif
    bool abort_on_out_of_gpu_memory = false;
}

const std::size_t Arena::align_size;

bool
Arena::isDeviceAccessible () const
{
#ifdef AMREX_USE_GPU
    return ! arena_info.use_cpu_memory;
#else
    return false;
#endif
}

bool
Arena::isHostAccessible () const
{
#ifdef AMREX_USE_GPU
    return (arena_info.use_cpu_memory ||
            arena_info.device_use_hostalloc ||
            arena_info.device_use_managed_memory);
#else
    return true;
#endif
}

bool
Arena::isManaged () const
{
#ifdef AMREX_USE_GPU
    return (! arena_info.use_cpu_memory)
        && (! arena_info.device_use_hostalloc)
        &&    arena_info.device_use_managed_memory;
#else
    return false;
#endif
}

bool
Arena::isDevice () const
{
#ifdef AMREX_USE_GPU
    return (! arena_info.use_cpu_memory)
        && (! arena_info.device_use_hostalloc)
        && (! arena_info.device_use_managed_memory);
#else
    return false;
#endif
}

bool
Arena::isPinned () const
{
#ifdef AMREX_USE_GPU
    return (! arena_info.use_cpu_memory)
        &&    arena_info.device_use_hostalloc;
#else
    return false;
#endif
}

bool
Arena::hasFreeDeviceMemory (std::size_t)
{
    return true;
}

void
Arena::registerForProfiling (const std::string&)
{
    amrex::Abort("Profiling is not implemented for this type of Arena");
}

std::size_t
Arena::align (std::size_t s)
{
    return amrex::aligned_size(align_size, s);
}

void*
Arena::allocate_system (std::size_t nbytes) // NOLINT(readability-make-member-function-const)
{
    void * p;
#ifdef AMREX_USE_GPU
    if (arena_info.use_cpu_memory)
    {
        p = std::malloc(nbytes);
#ifndef _WIN32
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        if (p && (nbytes > 0) && arena_info.device_use_hostalloc) mlock(p, nbytes);
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#endif
    }
    else if (arena_info.device_use_hostalloc)
    {
        AMREX_HIP_OR_CUDA_OR_SYCL(
            AMREX_HIP_SAFE_CALL (hipHostMalloc(&p, nbytes, hipHostMallocMapped|hipHostMallocNonCoherent));,
            AMREX_CUDA_SAFE_CALL(cudaHostAlloc(&p, nbytes, cudaHostAllocMapped));,
            p = sycl::malloc_host(nbytes, Gpu::Device::syclContext()));
    }
    else
    {
        std::size_t free_mem_avail = Gpu::Device::freeMemAvailable();
        if (nbytes >= free_mem_avail) {
            free_mem_avail += freeUnused_protected(); // For CArena, mutex has already acquired
            if (abort_on_out_of_gpu_memory && nbytes >= free_mem_avail) {
                amrex::Abort("Out of gpu memory. Free: " + std::to_string(free_mem_avail)
                             + " Asked: " + std::to_string(nbytes));
            }
        }

        if (arena_info.device_use_managed_memory)
        {
            AMREX_HIP_OR_CUDA_OR_SYCL
                (AMREX_HIP_SAFE_CALL(hipMallocManaged(&p, nbytes));,
                 AMREX_CUDA_SAFE_CALL(cudaMallocManaged(&p, nbytes));,
                 p = sycl::malloc_shared(nbytes, Gpu::Device::syclDevice(), Gpu::Device::syclContext()));
#ifdef AMREX_USE_HIP
            // Otherwise atomiAdd won't work because we instruct the compiler to do unsafe atomics
            AMREX_HIP_SAFE_CALL(hipMemAdvise(p, nbytes, hipMemAdviseSetCoarseGrain,
                                             Gpu::Device::deviceId()));
#endif
            if (arena_info.device_set_readonly)
            {
                Gpu::Device::mem_advise_set_readonly(p, nbytes);
            }
            if (arena_info.device_set_preferred)
            {
                const int device = Gpu::Device::deviceId();
                Gpu::Device::mem_advise_set_preferred(p, nbytes, device);
            }
        }
        else
        {
            AMREX_HIP_OR_CUDA_OR_SYCL
                (AMREX_HIP_SAFE_CALL ( hipMalloc(&p, nbytes));,
                 AMREX_CUDA_SAFE_CALL(cudaMalloc(&p, nbytes));,
                 p = sycl::malloc_device(nbytes, Gpu::Device::syclDevice(), Gpu::Device::syclContext()));
        }
    }
#else
    p = std::malloc(nbytes);
#ifndef _WIN32
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    if (p && (nbytes > 0) && arena_info.device_use_hostalloc) mlock(p, nbytes);
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#endif
#endif
    if (p == nullptr) amrex::Abort("Sorry, malloc failed");
    return p;
}

void
Arena::deallocate_system (void* p, std::size_t nbytes) // NOLINT(readability-make-member-function-const)
{
#ifdef AMREX_USE_GPU
    if (arena_info.use_cpu_memory)
    {
        if (p && arena_info.device_use_hostalloc) AMREX_MUNLOCK(p, nbytes);
        std::free(p);
    }
    else if (arena_info.device_use_hostalloc)
    {
        AMREX_HIP_OR_CUDA_OR_SYCL
            (AMREX_HIP_SAFE_CALL ( hipHostFree(p));,
             AMREX_CUDA_SAFE_CALL(cudaFreeHost(p));,
             sycl::free(p,Gpu::Device::syclContext()));
    }
    else
    {
        AMREX_HIP_OR_CUDA_OR_SYCL
            (AMREX_HIP_SAFE_CALL ( hipFree(p));,
             AMREX_CUDA_SAFE_CALL(cudaFree(p));,
             sycl::free(p,Gpu::Device::syclContext()));
    }
#else
    if (p && arena_info.device_use_hostalloc) AMREX_MUNLOCK(p, nbytes);
    std::free(p);
#endif
}

namespace {

    class NullArena final
        : public Arena
    {
        void* alloc (std::size_t) override { return nullptr; }
        void free (void*) override {}
    };

    Arena* The_Null_Arena ()
    {
        static NullArena the_null_arena;
        return &the_null_arena;
    }

    Arena* The_BArena ()
    {
        static BArena the_barena;
        return &the_barena;
    }
}

void
Arena::Initialize ()
{
    if (initialized) return;
    initialized = true;

    // see reason on allowed reuse of the default CPU BArena in Arena::Finalize
    BL_ASSERT(the_arena == nullptr || the_arena == The_BArena());
    BL_ASSERT(the_async_arena == nullptr);
    BL_ASSERT(the_device_arena == nullptr || the_device_arena == The_BArena());
    BL_ASSERT(the_managed_arena == nullptr || the_managed_arena == The_BArena());
    BL_ASSERT(the_pinned_arena == nullptr);
    BL_ASSERT(the_cpu_arena == nullptr || the_cpu_arena == The_BArena());

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    the_arena_init_size = 1024L*1024L*1024L; // xxxxx SYCL: todo
#else
    the_arena_init_size = Gpu::Device::totalGlobalMem() / Gpu::Device::numDevicePartners() / 4L * 3L;
#endif

    the_pinned_arena_release_threshold = Gpu::Device::totalGlobalMem();
#endif

    ParmParse pp("amrex");
    pp.queryAdd(        "the_arena_init_size",         the_arena_init_size);
    pp.queryAdd( "the_device_arena_init_size",  the_device_arena_init_size);
    pp.queryAdd("the_managed_arena_init_size", the_managed_arena_init_size);
    pp.queryAdd( "the_pinned_arena_init_size",  the_pinned_arena_init_size);
    pp.queryAdd(       "the_arena_release_threshold" ,         the_arena_release_threshold);
    pp.queryAdd( "the_device_arena_release_threshold",  the_device_arena_release_threshold);
    pp.queryAdd("the_managed_arena_release_threshold", the_managed_arena_release_threshold);
    pp.queryAdd( "the_pinned_arena_release_threshold",  the_pinned_arena_release_threshold);
    pp.queryAdd(  "the_async_arena_release_threshold",   the_async_arena_release_threshold);
    pp.queryAdd("the_arena_is_managed", the_arena_is_managed);
    pp.queryAdd("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);

    {
#if defined(BL_COALESCE_FABS) || defined(AMREX_USE_GPU)
        ArenaInfo ai{};
        ai.SetReleaseThreshold(the_arena_release_threshold);
        if (the_arena_is_managed) {
            the_arena = new CArena(0, ai.SetPreferred());
#ifdef AMREX_USE_GPU
            the_arena->registerForProfiling("Managed Memory");
#else
            the_arena->registerForProfiling("Cpu Memory");
#endif
        } else {
            the_arena = new CArena(0, ai.SetDeviceMemory());
#ifdef AMREX_USE_GPU
            the_arena->registerForProfiling("Device Memory");
#else
            the_arena->registerForProfiling("Cpu Memory");
#endif
        }
#ifdef AMREX_USE_GPU
        BL_PROFILE("The_Arena::Initialize()");
        void *p = the_arena->alloc(static_cast<std::size_t>(the_arena_init_size));
        the_arena->free(p);
#endif
#else
        the_arena = The_BArena();
#endif
    }

    the_async_arena = new PArena(the_async_arena_release_threshold);

#ifdef AMREX_USE_GPU
    if (the_arena->isDevice()) {
        the_device_arena = the_arena;
    } else {
        the_device_arena = new CArena(0, ArenaInfo{}.SetDeviceMemory().SetReleaseThreshold
                                      (the_device_arena_release_threshold));
        the_device_arena->registerForProfiling("Device Memory");
    }
#else
    the_device_arena = The_BArena();
#endif

#ifdef AMREX_USE_GPU
    if (the_arena->isManaged()) {
        the_managed_arena = the_arena;
    } else {
        the_managed_arena = new CArena(0, ArenaInfo{}.SetReleaseThreshold
                                       (the_managed_arena_release_threshold));
        the_managed_arena->registerForProfiling("Managed Memory");
    }
#else
    the_managed_arena = The_BArena();
#endif

    // When USE_CUDA=FALSE, we call mlock to pin the cpu memory.
    // When USE_CUDA=TRUE, we call cudaHostAlloc to pin the host memory.
    the_pinned_arena = new CArena(0, ArenaInfo{}.SetHostAlloc().SetReleaseThreshold
                                  (the_pinned_arena_release_threshold));
    the_pinned_arena->registerForProfiling("Pinned Memory");

    if (the_device_arena_init_size > 0 && the_device_arena != the_arena) {
        BL_PROFILE("The_Device_Arena::Initialize()");
        void *p = the_device_arena->alloc(the_device_arena_init_size);
        the_device_arena->free(p);
    }

    if (the_managed_arena_init_size > 0 && the_managed_arena != the_arena) {
        BL_PROFILE("The_Managed_Arena::Initialize()");
        void *p = the_managed_arena->alloc(the_managed_arena_init_size);
        the_managed_arena->free(p);
    }

    if (the_pinned_arena_init_size > 0) {
        BL_PROFILE("The_Pinned_Arena::Initialize()");
        void *p = the_pinned_arena->alloc(the_pinned_arena_init_size);
        the_pinned_arena->free(p);
    }

    the_cpu_arena = The_BArena();

    // Initialize the null arena
    auto* null_arena = The_Null_Arena();
    amrex::ignore_unused(null_arena);
}

void
Arena::PrintUsage ()
{
#ifdef AMREX_USE_GPU
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    {
        Long min_megabytes = Gpu::Device::totalGlobalMem() / (1024*1024);
        Long max_megabytes = min_megabytes;
        ParallelDescriptor::ReduceLongMin(min_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_megabytes, IOProc);
#ifdef AMREX_USE_MPI
        amrex::Print() << "Total GPU global memory (MB) spread across MPI: ["
                       << min_megabytes << " ... " << max_megabytes << "]\n";
#else
        amrex::Print() << "Total GPU global memory (MB): " << min_megabytes << "\n";
#endif
    }
    {
        Long min_megabytes = Gpu::Device::freeMemAvailable() / (1024*1024);
        Long max_megabytes = min_megabytes;
        ParallelDescriptor::ReduceLongMin(min_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_megabytes, IOProc);
#ifdef AMREX_USE_MPI
        amrex::Print() << "Free  GPU global memory (MB) spread across MPI: ["
                       << min_megabytes << " ... " << max_megabytes << "]\n";
#else
        amrex::Print() << "Free  GPU global memory (MB): " << min_megabytes << "\n";
#endif
    }
#endif
    if (The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Arena());
        if (p) {
            p->PrintUsage("The         Arena");
        }
    }
    if (The_Device_Arena() && The_Device_Arena() != The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Device_Arena());
        if (p) {
            p->PrintUsage("The  Device Arena");
        }
    }
    if (The_Managed_Arena() && The_Managed_Arena() != The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Managed_Arena());
        if (p) {
            p->PrintUsage("The Managed Arena");
        }
    }
    if (The_Pinned_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Pinned_Arena());
        if (p) {
            p->PrintUsage("The  Pinned Arena");
        }
    }
}

void
Arena::PrintUsageToFiles (const std::string& filename, const std::string& message)
{
    std::ofstream ofs(filename+"."+std::to_string(ParallelDescriptor::MyProc()),
                      std::ios_base::app);
    if (!ofs.is_open()) {
        amrex::Error("Could not open file for appending in amrex::Arena::PrintUsageToFiles()");
    }

    ofs << message << "\n";

#ifdef AMREX_USE_GPU
    Long megabytes = Gpu::Device::totalGlobalMem() / (1024*1024);
    ofs << "    Total GPU global memory (MB): " << megabytes << "\n";

    megabytes = Gpu::Device::freeMemAvailable() / (1024*1024);
    ofs << "    Free  GPU global memory (MB): " << megabytes << "\n";
#endif

    if (The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Arena());
        if (p) {
            p->PrintUsage(ofs, "The         Arena", "    ");
        }
    }
    if (The_Device_Arena() && The_Device_Arena() != The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Device_Arena());
        if (p) {
            p->PrintUsage(ofs, "The  Device Arena", "    ");
        }
    }
    if (The_Managed_Arena() && The_Managed_Arena() != The_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Managed_Arena());
        if (p) {
            p->PrintUsage(ofs, "The Managed Arena", "    ");
        }
    }
    if (The_Pinned_Arena()) {
        auto* p = dynamic_cast<CArena*>(The_Pinned_Arena());
        if (p) {
            p->PrintUsage(ofs, "The  Pinned Arena", "    ");
        }
    }

    ofs << "\n";
}

void
Arena::Finalize ()
{
#ifdef AMREX_USE_GPU
    if (amrex::Verbose() > 0) {
#else
    if (amrex::Verbose() > 1) {
#endif
        PrintUsage();
    }

    initialized = false;

    // we reset Arenas unless they are the default "CPU malloc/free" BArena
    // this is because we want to allow users to free their UB objects
    // that they forgot to destruct after amrex::Finalize():
    //   amrex::Initialize(...);
    //   MultiFab mf(...);  // this should be scoped in { ... }
    //   amrex::Finalize();
    // mf cannot be used now, but it can at least be freed without a segfault
    if (!dynamic_cast<BArena*>(the_device_arena)) {
        if (the_device_arena != the_arena) {
            delete the_device_arena;
        }
        the_device_arena = nullptr;
    }

    if (!dynamic_cast<BArena*>(the_managed_arena)) {
        if (the_managed_arena != the_arena) {
            delete the_managed_arena;
        }
        the_managed_arena = nullptr;
    }

    if (!dynamic_cast<BArena*>(the_arena)) {
        delete the_arena;
        the_arena = nullptr;
    }

    delete the_async_arena;
    the_async_arena = nullptr;

    delete the_pinned_arena;
    the_pinned_arena = nullptr;

    if (!dynamic_cast<BArena*>(the_cpu_arena)) {
        delete the_cpu_arena;
        the_cpu_arena = nullptr;
    }
}

Arena*
The_Arena ()
{
    if        (the_arena) {
        return the_arena;
    } else {
        return The_Null_Arena();
    }
}

Arena*
The_Async_Arena ()
{
    if        (the_async_arena) {
        return the_async_arena;
    } else {
        return The_Null_Arena();
    }
}

Arena*
The_Device_Arena ()
{
    if        (the_device_arena) {
        return the_device_arena;
    } else {
        return The_Null_Arena();
    }
}

Arena*
The_Managed_Arena ()
{
    if        (the_managed_arena) {
        return the_managed_arena;
    } else {
        return The_Null_Arena();
    }
}

Arena*
The_Pinned_Arena ()
{
    if        (the_pinned_arena) {
        return the_pinned_arena;
    } else {
        return The_Null_Arena();
    }
}

Arena*
The_Cpu_Arena ()
{
    if        (the_cpu_arena) {
        return the_cpu_arena;
    } else {
        return The_Null_Arena();
    }
}

}
