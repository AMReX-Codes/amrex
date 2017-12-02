#ifndef AMREX_INTCONV_H_
#define AMREX_INTCONV_H_

#include <iostream>
#include <cstring>
#include <stdint.h>

#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"

/// These are the integer types we support right now

template<typename T> T swapBytes(T value);

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

template <typename T>
void writeIntData(const int* data, std::size_t size, std::ostream& os,
                  const amrex::IntDescriptor& id)
{
    T value;
    bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
    for (std::size_t j = 0; j < size; ++j) {
        value = data[j];
        if (swapEndian) value = swapBytes<T>(value);
        os.write((char*) &value, sizeof(T));
    }
}

template <typename T>
void readIntData(int* data, std::size_t size, std::istream& is,
                 const amrex::IntDescriptor& id)
{
    T value;
    bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
    for (std::size_t j = 0; j < size; ++j) {
        is.read((char*) &value, sizeof(T));
        if (swapEndian) value = swapBytes<T>(value);
        data[j] = value;
    }
}

#endif

