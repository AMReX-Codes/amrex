
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

    amrex::ExecOnFinalize(Arena::Finalize);
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
        if (The_Cuda_Arena()) {
            CArena* p = dynamic_cast<CArena*>(The_Cuda_Arena());
            if (p) {
                long min_kilobytes = p->heap_space_used() / 1024;
                long max_kilobytes = min_kilobytes;
                ParallelDescriptor::ReduceLongMin(min_kilobytes, IOProc);
                ParallelDescriptor::ReduceLongMax(max_kilobytes, IOProc);
                amrex::Print() << "[The    CUDA Arena] space (kilobyte) used spread across MPI: ["
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
#endif
}

}
