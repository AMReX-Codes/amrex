//BL_COPYRIGHT_NOTICE

//
// $Id: BArena.cpp,v 1.6 2000-04-24 17:52:32 car Exp $
//

#include <BArena.H>
#include <BoxLib.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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

#ifdef BL_NAMESPACE
}
#endif

