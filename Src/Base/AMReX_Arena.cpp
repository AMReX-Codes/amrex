
#include <AMReX_Arena.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

#ifndef AMREX_FORTRAN_BOXLIB
#include <AMReX.H>
#endif

namespace amrex {

namespace {
    // The_Arena and The_Nvar_Arena are initialized in BaseFab
    Arena* the_cuda_arena = nullptr;
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
    delete the_cuda_arena;
    the_cuda_arena = nullptr;
    initialized = false;
#endif
}

Arena*
The_Cuda_Arena ()
{
    if (the_cuda_arena == nullptr) {
#if AMREX_USE_CUDA
        the_cuda_arena = new CArena;
        the_cuda_arena->SetDeviceMemory();
#else
        the_cuda_arena = new BArena;
#endif
    }

    return the_cuda_arena;
}

}
