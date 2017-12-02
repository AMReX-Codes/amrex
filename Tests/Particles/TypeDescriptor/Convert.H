#ifndef CONVERT_H_
#define CONVERT_H_

#include <iostream>
#include <cstring>
#include <stdint.h>

#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"

union Int32 {
    char bytes[4];
    std::int32_t n;
};

union Int64 {
    char bytes[8];
    std::int64_t n;
};

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
    return (val << 32) | (val >> 32);
}

template <typename T>
void writeIntData(const int* data, std::size_t size, std::ostream& os,
                  const amrex::IntDescriptor& id)
{
    T value;
    bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
    for (std::size_t j = 0; j < size; ++j) {
        value.n = data[j];
        // swap here;
        for (int i = 0; i < id.numBytes(); ++i) {
            os << value.bytes[i];
        }
    }
}

template <typename T>
void readIntData(int* data, std::size_t size, std::istream& is,
                 const amrex::IntDescriptor& id)
{
    T value;
    bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
    for (std::size_t j = 0; j < size; ++j) {
        for (int i = 0; i < id.numBytes(); ++i) {
            is.read((char*) &value.bytes[i], sizeof(char));
        }
        // swap here;
        data[j] = value.n;
    }
}

#endif

