//
// $Id: Arena.cpp,v 1.2 2001-07-19 16:57:30 lijewski Exp $
//

#include <Arena.H>
#include <BoxLib.H>

Arena::~Arena () {}

size_t
Arena::align (size_t s)
{
    size_t x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}
