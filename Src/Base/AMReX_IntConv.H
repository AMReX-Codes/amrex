#ifndef AMREX_INTCONV_H_
#define AMREX_INTCONV_H_

#include <iostream>
#include <cstring>
#include <stdint.h>

#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"

namespace amrex {

    template<typename T> T swapBytes(T value);
    
    template <typename To, typename From>
    void writeIntData(const From* data, std::size_t size, std::ostream& os,
                      const amrex::IntDescriptor& id)
    {
        To value;
        bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
        for (std::size_t j = 0; j < size; ++j) {
            value = data[j];
            if (swapEndian) value = swapBytes<To>(value);
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
            if (swapEndian) value = swapBytes<From>(value);
            data[j] = value;
        }
    }   
}

#endif

