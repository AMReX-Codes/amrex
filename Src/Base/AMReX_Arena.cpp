
#include <AMReX_Arena.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>
#include <AMReX_DArena.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>

namespace amrex {

namespace {
    bool initialized = false;

    Arena* the_arena = nullptr;
    Arena* the_device_arena = nullptr;
    Arena* the_managed_arena = nullptr;
    Arena* the_pinned_arena = nullptr;

    bool use_buddy_allocator = false;
    long buddy_allocator_size = 0L;
}

const unsigned int Arena::align_size;

Arena::~Arena () {}

std::size_t
Arena::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

void*
Arena::allocate_system (std::size_t nbytes)
{
#ifdef AMREX_USE_CUDA
    void * p;
    if (arena_info.device_use_hostalloc)
    {
        AMREX_GPU_SAFE_CALL(cudaHostAlloc(&p, nbytes, cudaHostAllocMapped));
    }
    else if (arena_info.device_use_managed_memory)
    {
        AMREX_GPU_SAFE_CALL(cudaMallocManaged(&p, nbytes));
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
        AMREX_GPU_SAFE_CALL(cudaMalloc(&p, nbytes));
    }
    return p;
#else
    return std::malloc(nbytes);
#endif
}

void
Arena::deallocate_system (void* p)
{
#ifdef AMREX_USE_CUDA
    if (arena_info.device_use_hostalloc)
    {
        AMREX_GPU_SAFE_CALL(cudaFreeHost(p));
    }
    else
    {
        AMREX_GPU_SAFE_CALL(cudaFree(p));
    }
#else
    std::free(p);
#endif
}

void
Arena::Initialize ()
{
    if (initialized) return;
    initialized = true;

    BL_ASSERT(the_arena == nullptr);
    BL_ASSERT(the_device_arena == nullptr);
    BL_ASSERT(the_managed_arena == nullptr);
    BL_ASSERT(the_pinned_arena == nullptr);

    ParmParse pp("amrex");
    pp.query("use_buddy_allocator", use_buddy_allocator);
    pp.query("buddy_allocator_size", buddy_allocator_size);

#ifdef AMREX_USE_GPU
    if (use_buddy_allocator)
    {
        if (buddy_allocator_size <= 0) {
            buddy_allocator_size = Gpu::Device::totalGlobalMem() / 4 * 3;
        }
        std::size_t chunk = 512*1024*1024;
        buddy_allocator_size = (buddy_allocator_size/chunk) * chunk;
        the_arena = new DArena(buddy_allocator_size, 512, ArenaInfo().SetPreferred());
    }
    else
#endif
    {
#if defined(BL_COALESCE_FABS) || defined(AMREX_USE_GPU)
        the_arena = new CArena(0, ArenaInfo().SetPreferred());
#else
        the_arena = new BArena;
#endif
    }

#ifdef AMREX_USE_GPU
    the_device_arena = new CArena(0, ArenaInfo().SetDeviceMemory());
#else
    the_device_arena = new BArena;
#endif

#if defined(AMREX_USE_GPU)
    the_managed_arena = new CArena;
#else
    the_managed_arena = new BArena;
#endif

#if defined(AMREX_USE_GPU)
//    const std::size_t hunk_size = 64 * 1024;
//    the_pinned_arena = new CArena(hunk_size);
    the_pinned_arena = new CArena(0, ArenaInfo().SetHostAlloc());
#else
    the_pinned_arena = new BArena;
#endif

    std::size_t N = 1024*1024*8;
    void *p = the_arena->alloc(N);
    the_arena->free(p);

    p = the_device_arena->alloc(N);
    the_device_arena->free(p);

    p = the_managed_arena->alloc(N);
    the_managed_arena->free(p);

    p = the_pinned_arena->alloc(N);
    the_pinned_arena->free(p);
}

void
Arena::PrintUsage ()
{
    if (amrex::Verbose() > 0) {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        if (The_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
#ifdef AMREX_USE_MPI
                amrex::Print() << "[The         Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
#else
                amrex::Print() << "[The         Arena] space (kilobyte): " << min_kilobytes << "\n";
#endif
            }
        }
        if (The_Device_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Device_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
#ifdef AMREX_USE_MPI
                amrex::Print() << "[The  Device Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
#else
                amrex::Print() << "[The  Device Arena] space (kilobyte): " << min_kilobytes << "\n";
#endif
            }
        }
        if (The_Managed_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Managed_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
#ifdef AMREX_USE_MPI
                amrex::Print() << "[The Managed Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
#else
                amrex::Print() << "[The Managed Arena] space (kilobyte): " << min_kilobytes << "\n";
#endif
            }
        }
        if (The_Pinned_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Pinned_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
#ifdef AMREX_USE_MPI
                amrex::Print() << "[The  Pinned Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
#else
                amrex::Print() << "[The  Pinned Arena] space (kilobyte): " << min_kilobytes << "\n";
#endif
            }
        }
    }
}
    
void
Arena::Finalize ()
{
    PrintUsage();
    
    initialized = false;
    
    delete the_arena;
    the_arena = nullptr;
    
    delete the_device_arena;
    the_device_arena = nullptr;
    
    delete the_managed_arena;
    the_managed_arena = nullptr;
    
    delete the_pinned_arena;
    the_pinned_arena = nullptr;
}
    
Arena*
The_Arena ()
{
    BL_ASSERT(the_arena != nullptr);
    return the_arena;
}

Arena*
The_Device_Arena ()
{
    BL_ASSERT(the_device_arena != nullptr);
    return the_device_arena;
}

Arena*
The_Managed_Arena ()
{
    BL_ASSERT(the_managed_arena != nullptr);
    return the_managed_arena;
}

Arena*
The_Pinned_Arena ()
{
    BL_ASSERT(the_pinned_arena != nullptr);
    return the_pinned_arena;
}

}
