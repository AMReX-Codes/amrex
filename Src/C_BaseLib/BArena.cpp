//BL_COPYRIGHT_NOTICE

//
// $Id: BArena.cpp,v 1.5 1998-02-09 20:46:52 lijewski Exp $
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
