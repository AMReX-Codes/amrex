//BL_COPYRIGHT_NOTICE

//
// $Id: BArena.cpp,v 1.4 1998-02-06 18:22:03 lijewski Exp $
//

#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (size_t _sz,
               void** _pt)
{
    void* pt = ::operator new(_sz);
    if (_pt != 0)
        *_pt = pt;
    return pt;
}

void BArena::free (void* pt)
{
    ::operator delete(pt);
}
