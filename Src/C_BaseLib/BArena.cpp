
#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (std::size_t _sz)
{
    return ::operator new(_sz);
}

void
BArena::free (void* pt)
{
    ::operator delete(pt);
}
