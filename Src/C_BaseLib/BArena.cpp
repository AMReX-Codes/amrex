//BL_COPYRIGHT_NOTICE

//
// $Id: BArena.cpp,v 1.2 1997-09-15 23:33:19 lijewski Exp $
//

#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (size_t _sz,
               void** _pt)
{
    void* pt = new char[_sz];
    if (pt == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    if (_pt != 0)
        *_pt = pt;
    return pt;
}

void BArena::free (void* pt)
{
    delete [] pt;
}
