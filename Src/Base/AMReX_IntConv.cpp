#include "AMReX_IntConv.H"

namespace amrex {

template<> std::int16_t swapBytes<std::int16_t>(std::int16_t val)
{
    return (val << 8) | ((val >> 8) & 0xFF);
}

template<> std::int32_t swapBytes<std::int32_t>(std::int32_t val)
{
    val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
    return (val << 16) | ((val >> 16) & 0xFFFF);
}

template<> std::int64_t swapBytes<std::int64_t>(std::int64_t val)
{
    val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
    val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
    return (val << 32) | ((val >> 32) & 0xFFFFFFFFULL);
}

template<> std::uint16_t swapBytes<std::uint16_t>(std::uint16_t val)
{
    return (val << 8) | (val >> 8 );
}

template<> std::uint32_t swapBytes<std::uint32_t>(std::uint32_t val)
{
    val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF ); 
    return (val << 16) | (val >> 16);
}

template<> std::uint64_t swapBytes<std::uint64_t>(std::uint64_t val)
{
    val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
    val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
    return (val << 32) | (val >> 32);
}

}

