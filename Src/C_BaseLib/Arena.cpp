//
// $Id: Arena.cpp,v 1.3 2004-01-07 21:18:19 car Exp $
//

#include <Arena.H>
#include <BoxLib.H>

Arena::~Arena () {}

std::size_t
Arena::align (std::size_t s)
{
    std::size_t x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}
