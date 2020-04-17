#ifndef AMREX_VECTORIO_H_
#define AMREX_VECTORIO_H_

#include <iostream>

#include "AMReX_FPC.H"
#include "AMReX_FabConv.H"
#include "AMReX_IntConv.H"

namespace amrex {

/**
* \brief Functions for writing integer data to disk in a portable, self-describing manner.
*/

    //! Write int data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! IntDescriptor that describes the data format to use for writing. If no
    //! IntDescriptor is provided, the data will be written using the native
    //! format for your machine.
    void writeIntData(const int* data, std::size_t size, std::ostream& os,
                      const IntDescriptor& id = FPC::NativeIntDescriptor());

    //! Read int data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and an IntDescriptor
    //! that describes the format of the data on disk. 
    //! The buffer is assumed to be large enough to store 'size' integers, and it is
    //! the user's reponsiblity to allocate this data.
    void readIntData(int* data, std::size_t size, std::istream& is, const IntDescriptor& id);

    //! Write long data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! IntDescriptor that describes the data format to use for writing. If no
    //! IntDescriptor is provided, the data will be written using the native
    //! format for your machine.
    void writeLongData(const Long* data, std::size_t size, std::ostream& os,
                       const IntDescriptor& id = FPC::NativeLongDescriptor());

    //! Read int data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and an IntDescriptor
    //! that describes the format of the data on disk. 
    //! The buffer is assumed to be large enough to store 'size' longs, and it is
    //! the user's reponsiblity to allocate this data.
    void readLongData(Long* data, std::size_t size, std::istream& is, const IntDescriptor& id);

    //! Write Real data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.    
    void writeRealData(const Real* data, std::size_t size, std::ostream& os,
                       const RealDescriptor& rd = FPC::NativeRealDescriptor());

    //! Read Real data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    void readRealData(Real* data, std::size_t size, std::istream& is,
                      const RealDescriptor& rd);

    //! Write float data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.
    void writeFloatData(const float* data, std::size_t size, std::ostream& os,
                        const RealDescriptor& rd = FPC::Native32RealDescriptor());

    //! Read float data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    void readFloatData(float* data, std::size_t size, std::istream& is,
                       const RealDescriptor& rd);

    //! Write double data to the ostream. The arguments are a pointer to data
    //! to write, the size of the data buffer, the ostream, and an optional
    //! RealDescriptor that describes the data format to use for writing. If no
    //! RealDescriptor is provided, the data will be written using the native
    //! format for your machine.
    void writeDoubleData(const double* data, std::size_t size, std::ostream& os,
                         const RealDescriptor& rd = FPC::Native64RealDescriptor());

    //! Read double data from the istream. The arguments are a pointer to data buffer
    //! to read into, the size of that buffer, the istream, and a RealDescriptor
    //! that describes the format of the data on disk.
    //! The buffer is assumed to be large enough to store 'size' Reals, and it is
    //! the user's reponsiblity to allocate this data.
    void readDoubleData(double* data, std::size_t size, std::istream& is,
                        const RealDescriptor& rd);

    inline void writeData (int const* data, std::size_t size, std::ostream& os)
    {
        writeIntData(data, size, os);
    }

    inline void writeData (Long const* data, std::size_t size, std::ostream& os)
    {
        writeLongData(data, size, os);
    }

    inline void writeData (float const* data, std::size_t size, std::ostream& os)
    {
        writeFloatData(data, size, os);
    }

    inline void writeData (double const* data, std::size_t size, std::ostream& os)
    {
        writeDoubleData(data, size, os);
    }

    inline void readData (int * data, std::size_t size, std::istream& is)
    {
        readIntData(data, size, is, FPC::NativeIntDescriptor());
    }

    inline void readData (Long * data, std::size_t size, std::istream& is)
    {
        readLongData(data, size, is, FPC::NativeLongDescriptor());
    }

    inline void readData (float * data, std::size_t size, std::istream& is)
    {
        readFloatData(data, size, is, FPC::Native32RealDescriptor());
    }

    inline void readData (double * data, std::size_t size, std::istream& is)
    {
        readDoubleData(data, size, is, FPC::Native64RealDescriptor());
    }
}

#endif 
