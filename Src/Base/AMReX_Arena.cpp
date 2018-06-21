
#include <AMReX_Arena.H>
#include <AMReX.H>

namespace amrex {

namespace {
#ifdef BL_COALESCE_FABS
  static std::unique_ptr<CArena> the_arena;
#else
  static std::unique_ptr<BArena> the_arena;
#endif
}

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

#ifdef BL_COALESCE_FABS
CArena*
#else
BArena*
#endif
The_Arena ()
{
    BL_ASSERT(the_arena != 0);

    return the_arena;
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
