
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
