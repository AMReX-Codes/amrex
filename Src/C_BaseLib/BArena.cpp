//
// $Id: BArena.cpp,v 1.10 2001-07-31 17:56:24 lijewski Exp $
//
#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (size_t _sz)
{
    return ::operator new(_sz);
}

void
BArena::free (void* pt)
{
    ::operator delete(pt);
}
