
#include <memory>
#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

namespace amrex {

namespace {
  static std::unique_ptr<Arena> the_arena;
}

#if !defined(AMREX_FORTRAN_BOXLIB)
void
Arena_Initialize()
{
#ifdef BL_COALESCE_FABS
    the_arena.reset(new CArena);
#else
    the_arena.reset(new BArena);
#endif

    amrex::ExecOnFinalize(amrex::Arena_Finalize);
}

void
Arena_Finalize()
{
    the_arena.reset();
}
#endif

Arena*
The_Arena ()
{
    BL_ASSERT(the_arena != 0);

    return the_arena.get();
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

}
