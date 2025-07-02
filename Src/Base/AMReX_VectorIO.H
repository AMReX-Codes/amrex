#ifndef AMREX_VECTORIO_H_
#define AMREX_VECTORIO_H_
#include <AMReX_Config.H>

#include <AMReX_FPC.H>
#include <AMReX_FabConv.H>
#include <AMReX_IntConv.H>

#include <iostream>

namespace amrex {

/**
* \brief Functions for writing integer data to disk in a portable, self-describing manner.
*/

    //! Write int data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! IntDescriptor that describes the data format to use for writing. If no
    //! IntDescriptor is provided, the data will be written using the native
    //! format for your machine.
    template <typename OS>
    void writeIntData (const int* data, std::size_t size, OS& os,
                       const IntDescriptor& id = FPC::NativeIntDescriptor())
    {
        if (id == FPC::NativeIntDescriptor())
        {
            os.write((char*) data, static_cast<std::streamsize>(size*sizeof(int)));
        }
        else if (id.numBytes() == 2)
        {
            writeIntData<std::int16_t, int>(data, size, os, id);
        }
        else if (id.numBytes() == 4)
        {
            writeIntData<std::int32_t, int>(data, size, os, id);
        }
        else if (id.numBytes() == 8)
        {
            writeIntData<std::int64_t, int>(data, size, os, id);
        }
        else {
            amrex::Error("Don't know how to work with this integer type.");
        }
    }

    //! Read int data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and an IntDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' integers, and it is
    //! the user's reponsiblity to allocate this data.
    template <typename IS>
    void readIntData (int* data, std::size_t size, IS& is, const IntDescriptor& id)
    {
        if (id == FPC::NativeIntDescriptor())
        {
            is.read((char*) data, static_cast<std::streamsize>(size * id.numBytes()));
        }
        else if (id.numBytes() == 2)
        {
            readIntData<int, std::int16_t>(data, size, is, id);
        }
        else if (id.numBytes() == 4)
        {
            readIntData<int, std::int32_t>(data, size, is, id);
        }
        else if (id.numBytes() == 8)
        {
            readIntData<int, std::int64_t>(data, size, is, id);
        }
        else {
            amrex::Error("Don't know how to work with this integer type.");
        }
    }

    //! Write long data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! IntDescriptor that describes the data format to use for writing. If no
    //! IntDescriptor is provided, the data will be written using the native
    //! format for your machine.
    template <typename OS>
    void writeLongData (const Long* data, std::size_t size, OS& os,
                        const IntDescriptor& id = FPC::NativeLongDescriptor())
    {
        if (id == FPC::NativeLongDescriptor())
        {
            os.write((char*) data, static_cast<std::streamsize>(size*sizeof(Long)));
        }
        else if (id.numBytes() == 2)
        {
            writeIntData<std::int16_t, Long>(data, size, os, id);
        }
        else if (id.numBytes() == 4)
        {
            writeIntData<std::int32_t, Long>(data, size, os, id);
        }
        else if (id.numBytes() == 8)
        {
            writeIntData<std::int64_t, Long>(data, size, os, id);
        }
        else {
            amrex::Error("Don't know how to work with this long type.");
        }
    }

    //! Read int data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and an IntDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' longs, and it is
    //! the user's reponsiblity to allocate this data.
    template <typename IS>
    void readLongData (Long* data, std::size_t size, IS& is, const IntDescriptor& id)
    {
        if (id == FPC::NativeLongDescriptor())
        {
            is.read((char*) data, static_cast<std::streamsize>(size * id.numBytes()));
        }
        else if (id.numBytes() == 2)
        {
            readIntData<Long, std::int16_t>(data, size, is, id);
        }
        else if (id.numBytes() == 4)
        {
            readIntData<Long, std::int32_t>(data, size, is, id);
        }
        else if (id.numBytes() == 8)
        {
            readIntData<Long, std::int64_t>(data, size, is, id);
        }
        else {
            amrex::Error("Don't know how to work with this long type.");
        }
    }

    //! Write Real data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.
    template <typename OS>
    void writeRealData (const Real* data, std::size_t size, OS& os,
                        const RealDescriptor& rd = FPC::NativeRealDescriptor())
    {
        RealDescriptor::convertFromNativeFormat(os, static_cast<Long>(size), data, rd);
    }

    //! Read Real data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    template <typename IS>
    void readRealData (Real* data, std::size_t size, IS& is,
                       const RealDescriptor& rd)
    {
        RealDescriptor::convertToNativeFormat(data, static_cast<Long>(size), is, rd);
    }

    //! Write float data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.
    template <typename OS>
    void writeFloatData (const float* data, std::size_t size, OS& os,
                         const RealDescriptor& rd = FPC::Native32RealDescriptor())
    {
        RealDescriptor::convertFromNativeFloatFormat(os, static_cast<Long>(size), data, rd);
    }

    //! Read float data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    template <typename IS>
    void readFloatData(float* data, std::size_t size, IS& is,
                       const RealDescriptor& rd)
    {
        RealDescriptor::convertToNativeFloatFormat(data, static_cast<Long>(size), is, rd);
    }

    //! Write double data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.
    template <typename OS>
    void writeDoubleData (const double* data, std::size_t size, OS& os,
                          const RealDescriptor& rd = FPC::Native64RealDescriptor())
    {
        RealDescriptor::convertFromNativeDoubleFormat(os, static_cast<Long>(size), data, rd);
    }

    //! Read double data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    template <typename IS>
    void readDoubleData (double* data, std::size_t size, IS& is,
                         const RealDescriptor& rd)
    {
        RealDescriptor::convertToNativeDoubleFormat(data, static_cast<Long>(size), is, rd);
    }

    template <typename OS>
    void writeData (int const* data, std::size_t size, OS& os)
    {
        writeIntData(data, size, os);
    }

    template <typename OS>
    void writeData (Long const* data, std::size_t size, OS& os)
    {
        writeLongData(data, size, os);
    }

    template <typename OS>
    void writeData (float const* data, std::size_t size, OS& os)
    {
        writeFloatData(data, size, os);
    }

    template <typename OS>
    void writeData (double const* data, std::size_t size, OS& os)
    {
        writeDoubleData(data, size, os);
    }

    template <typename IS>
    void readData (int * data, std::size_t size, IS& is)
    {
        readIntData(data, size, is, FPC::NativeIntDescriptor());
    }

    template <typename IS>
    void readData (Long * data, std::size_t size, IS& is)
    {
        readLongData(data, size, is, FPC::NativeLongDescriptor());
    }

    template <typename IS>
    void readData (float * data, std::size_t size, IS& is)
    {
        readFloatData(data, size, is, FPC::Native32RealDescriptor());
    }

    template <typename IS>
    void readData (double * data, std::size_t size, IS& is)
    {
        readDoubleData(data, size, is, FPC::Native64RealDescriptor());
    }
}

#endif
