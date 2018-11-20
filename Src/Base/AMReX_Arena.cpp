
#include <AMReX_Arena.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

#ifndef AMREX_FORTRAN_BOXLIB
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#endif

namespace amrex {

namespace {
    bool initialized = false;

    Arena* the_arena = nullptr;
    Arena* the_device_arena = nullptr;
    Arena* the_managed_arena = nullptr;
    Arena* the_pinned_arena = nullptr;
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

void
Arena::Initialize ()
{
#ifndef AMREX_FORTRAN_BOXLIB
    if (initialized) return;
    initialized = true;

    BL_ASSERT(the_arena == nullptr);
    BL_ASSERT(the_device_arena == nullptr);
    BL_ASSERT(the_managed_arena == nullptr);
    BL_ASSERT(the_pinned_arena == nullptr);
    
#if defined(BL_COALESCE_FABS)
    the_arena = new CArena;
#else
    the_arena = new BArena;
#endif
    
#ifdef AMREX_USE_GPU
    the_arena->SetPreferred();
#endif

#if AMREX_USE_GPU
    the_device_arena = new CArena;
    the_device_arena->SetDeviceMemory();
#else
    the_device_arena = new BArena;
#endif
    
#if defined(AMREX_USE_GPU)
    the_managed_arena = new CArena;
#else
    the_managed_arena = new BArena;
#endif

#if defined(AMREX_USE_GPU)
    const std::size_t hunk_size = 64 * 1024;
    the_pinned_arena = new CArena(hunk_size);
    the_pinned_arena->SetHostAlloc();
#else
    the_pinned_arena = new BArena;
#endif

#endif
}

void
Arena::Finalize ()
{
#ifndef AMREX_FORTRAN_BOXLIB
    initialized = false;
    if (amrex::Verbose() > 0) {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        if (The_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
                amrex::Print() << "[The         Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
            }
        }
        if (The_Device_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Device_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
                amrex::Print() << "[The  Device Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
            }
        }
        if (The_Managed_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Managed_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
                amrex::Print() << "[The Managed Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
            }
        }
        if (The_Pinned_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Pinned_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
                amrex::Print() << "[The  Pinned Arena] space (kilobyte) used spread across MPI: ["
                               << min_kilobytes << " ... " << max_kilobytes << "]\n";
            }
        }
    }

    delete the_arena;
    the_arena = nullptr;

    delete the_device_arena;
    the_device_arena = nullptr;

    delete the_managed_arena;
    the_managed_arena = nullptr;

    delete the_pinned_arena;
    the_pinned_arena = nullptr;
#endif
}

#ifndef AMREX_FORTRAN_BOXLIB

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

#endif
}
