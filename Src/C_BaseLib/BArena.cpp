
#include <BArena.H>
#include <BoxLib.H>
#include <ParallelDescriptor.H>

#ifdef BL_USE_UPCXX
#include <upcxx.h>
#endif

void*
BArena::alloc (std::size_t _sz)
{
#ifdef BL_USE_UPCXX
    return upcxx::allocate(_sz);
#else
    return ::operator new(_sz);
#endif
}

void
BArena::free (void* pt)
{
#ifdef BL_USE_UPCXX
    upcxx::deallocate(pt);
#else
    ::operator delete(pt);
#endif
}
