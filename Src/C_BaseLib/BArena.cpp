//BL_COPYRIGHT_NOTICE

//
// $Id: BArena.cpp,v 1.3 1997-12-11 23:25:39 lijewski Exp $
//

#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (size_t _sz,
               void** _pt)
{
    void* pt = new char[_sz];
    if (_pt != 0)
        *_pt = pt;
    return pt;
}

void BArena::free (void* pt)
{
    delete [] pt;
}
