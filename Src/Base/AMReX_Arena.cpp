
#include <AMReX_Arena.H>
#include <AMReX.H>

const unsigned int amrex::Arena::align_size;

amrex::Arena::~Arena () {}

std::size_t
amrex::Arena::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}
