
//
// $Id: BArena.cpp,v 1.7 2000-10-02 20:52:32 lijewski Exp $
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

