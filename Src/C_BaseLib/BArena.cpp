//
// $Id: BArena.cpp,v 1.11 2004-01-07 21:18:19 car Exp $
//
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
