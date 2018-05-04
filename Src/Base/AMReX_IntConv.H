#ifndef AMREX_INTCONV_H_
#define AMREX_INTCONV_H_

#include <iostream>
#include <cstring>
#include <cstdint>

#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"

namespace amrex {

    std::int16_t swapBytes(std::int16_t value);
    std::int32_t swapBytes(std::int32_t value);
    std::int64_t swapBytes(std::int64_t value);

    std::uint16_t swapBytes(std::uint16_t value);
    std::uint32_t swapBytes(std::uint32_t value);
    std::uint64_t swapBytes(std::uint64_t value);

    template <typename To, typename From>
    void writeIntData(const From* data, std::size_t size, std::ostream& os,
                      const amrex::IntDescriptor& id)
    {
        To value;
        bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
        for (std::size_t j = 0; j < size; ++j) {
            value = data[j];
            if (swapEndian) value = swapBytes(value);
            os.write((char*) &value, sizeof(To));
        }
    }
    
    template <typename To, typename From>
    void readIntData(To* data, std::size_t size, std::istream& is,
                     const amrex::IntDescriptor& id)
    {
        From value;
        bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
        for (std::size_t j = 0; j < size; ++j) {
            is.read((char*) &value, sizeof(From));
            if (swapEndian) value = swapBytes(value);
            data[j] = value;
        }
    }   
}

#endif

