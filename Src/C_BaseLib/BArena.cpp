//
// $Id: BArena.cpp,v 1.9 2001-07-19 16:57:30 lijewski Exp $
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
