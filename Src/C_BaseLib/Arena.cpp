//
// $Id: Arena.cpp,v 1.1 2001-07-19 15:26:47 lijewski Exp $
//

#include <Arena.H>
#include <BoxLib.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

Arena::~Arena () {}

size_t
Arena::align (size_t s)
{
    size_t x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}

#ifdef BL_NAMESPACE
}
#endif

