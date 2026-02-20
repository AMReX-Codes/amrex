#ifndef AMREX_INTCONV_H_
#define AMREX_INTCONV_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_IntConv.H
 * \brief Integer byte-swap and endian-aware I/O helpers.
 */

#include <AMReX_FPC.H>
#include <AMReX_FabConv.H>

#include <iostream>
#include <cstring>
#include <cstdint>

namespace amrex {

    //! Reverse the byte order of a 16-bit signed integer.
    std::int16_t swapBytes (std::int16_t value);
    //! Reverse the byte order of a 32-bit signed integer.
    std::int32_t swapBytes (std::int32_t value);
    //! Reverse the byte order of a 64-bit signed integer.
    std::int64_t swapBytes (std::int64_t value);

    //! Reverse the byte order of a 16-bit unsigned integer.
    std::uint16_t swapBytes (std::uint16_t value);
    //! Reverse the byte order of a 32-bit unsigned integer.
    std::uint32_t swapBytes (std::uint32_t value);
    //! Reverse the byte order of a 64-bit unsigned integer.
    std::uint64_t swapBytes (std::uint64_t value);

    /**
     * \brief Write \p size integer values from \p data to \p os,
     *        converting from type \p From to type \p To and
     *        byte-swapping if the endianness of \p id differs from native.
     */
    template <typename To, typename From>
    void writeIntData (const From* data, std::size_t size, std::ostream& os,
                       const amrex::IntDescriptor& id)
    {
        To value;
        bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
        for (std::size_t j = 0; j < size; ++j) {
            value = static_cast<To>(data[j]);
            if (swapEndian) { value = swapBytes(value); }
            os.write((char*) &value, sizeof(To));
        }
    }

    /**
     * \brief Read \p size integer values from \p is into \p data,
     *        converting from type \p From to type \p To and
     *        byte-swapping if the endianness of \p id differs from native.
     */
    template <typename To, typename From>
    void readIntData (To* data, std::size_t size, std::istream& is,
                      const amrex::IntDescriptor& id)
    {
        From value;
        bool swapEndian = (id.order() != amrex::FPC::NativeIntDescriptor().order());
        for (std::size_t j = 0; j < size; ++j) {
            is.read((char*) &value, sizeof(From));
            if (swapEndian) { value = swapBytes(value); }
            data[j] = static_cast<To>(value);
        }
    }
}

#endif
