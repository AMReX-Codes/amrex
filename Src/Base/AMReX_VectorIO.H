#ifndef AMREX_VECTORIO_H_
#define AMREX_VECTORIO_H_
#include <AMReX_Config.H>

#include <AMReX_FPC.H>
#include <AMReX_FabConv.H>
#include <AMReX_IntConv.H>

#include <iosfwd>

/**
 * \file AMReX_VectorIO.H
 * \brief Portable binary I/O helpers for int, long, and Real arrays.
 */

namespace amrex {

/**
 * \brief Functions for writing integer data to disk in a portable, self-describing manner.
 */

    /**
     * \brief Write int data to a stream using the provided descriptor.
     * \param data pointer to the beginning of the integer buffer.
     * \param size number of ints in the buffer.
     * \param os   destination stream.
     * \param id   descriptor describing the on-disk integer format (defaults to native).
     */
    void writeIntData (const int* data, std::size_t size, std::ostream& os,
                       const IntDescriptor& id = FPC::NativeIntDescriptor());

    /**
     * \brief Read int data from a stream using the provided descriptor.
     * \param data destination buffer (length \p size).
     * \param size number of ints to read.
     * \param is   source stream.
     * \param id   descriptor describing the on-disk integer format.
     */
    void readIntData (int* data, std::size_t size, std::istream& is, const IntDescriptor& id);

    /**
     * \brief Write long data to a stream using the provided descriptor.
     * \param data pointer to the beginning of the long buffer.
     * \param size number of longs in the buffer.
     * \param os   destination stream.
     * \param id   descriptor describing the on-disk integer format (defaults to native).
     */
    void writeLongData (const Long* data, std::size_t size, std::ostream& os,
                        const IntDescriptor& id = FPC::NativeLongDescriptor());

    /**
     * \brief Read long data from a stream using the provided descriptor.
     * \param data destination buffer (length \p size).
     * \param size number of longs to read.
     * \param is   source stream.
     * \param id   descriptor describing the on-disk integer format.
     */
    void readLongData (Long* data, std::size_t size, std::istream& is, const IntDescriptor& id);

    /**
     * \brief Write Real data to a stream using the provided descriptor.
     * \param data pointer to the beginning of the Real buffer.
     * \param size number of Reals in the buffer.
     * \param os   destination stream.
     * \param rd   descriptor describing the on-disk floating-point format (defaults to native).
     */
    void writeRealData (const Real* data, std::size_t size, std::ostream& os,
                        const RealDescriptor& rd = FPC::NativeRealDescriptor());

    /**
     * \brief Read Real data from a stream using the provided descriptor.
     * \param data destination buffer (length \p size).
     * \param size number of Reals to read.
     * \param is   source stream.
     * \param rd   descriptor describing the on-disk floating-point format.
     */
    void readRealData (Real* data, std::size_t size, std::istream& is,
                       const RealDescriptor& rd);

    /**
     * \brief Write float data to a stream using the provided descriptor.
     * \param data pointer to the beginning of the float buffer.
     * \param size number of floats in the buffer.
     * \param os   destination stream.
     * \param rd   descriptor describing the on-disk floating-point format (defaults to native).
     */
    void writeFloatData (const float* data, std::size_t size, std::ostream& os,
                         const RealDescriptor& rd = FPC::Native32RealDescriptor());

    /**
     * \brief Read float data from a stream using the provided descriptor.
     * \param data destination buffer (length \p size).
     * \param size number of floats to read.
     * \param is   source stream.
     * \param rd   descriptor describing the on-disk floating-point format.
     */
    void readFloatData(float* data, std::size_t size, std::istream& is,
                       const RealDescriptor& rd);

    /**
     * \brief Write double data to a stream using the provided descriptor.
     * \param data pointer to the beginning of the double buffer.
     * \param size number of doubles in the buffer.
     * \param os   destination stream.
     * \param rd   descriptor describing the on-disk floating-point format (defaults to native).
     */
    void writeDoubleData (const double* data, std::size_t size, std::ostream& os,
                          const RealDescriptor& rd = FPC::Native64RealDescriptor());

    /**
     * \brief Read double data from a stream using the provided descriptor.
     * \param data destination buffer (length \p size).
     * \param size number of doubles to read.
     * \param is   source stream.
     * \param rd   descriptor describing the on-disk floating-point format.
     */
    void readDoubleData (double* data, std::size_t size, std::istream& is,
                         const RealDescriptor& rd);

    /**
     * \brief Convenience overload forwarding to writeIntData() with the native descriptor.
     * \param data pointer to the integer buffer.
     * \param size number of ints to write.
     * \param os   destination stream.
     */
    inline void writeData (int const* data, std::size_t size, std::ostream& os)
    {
        writeIntData(data, size, os);
    }

    /**
     * \brief Convenience overload forwarding to writeLongData() with the native descriptor.
     * \param data pointer to the long buffer.
     * \param size number of longs to write.
     * \param os   destination stream.
     */
    inline void writeData (Long const* data, std::size_t size, std::ostream& os)
    {
        writeLongData(data, size, os);
    }

    /**
     * \brief Convenience overload forwarding to writeFloatData() with the native descriptor.
     * \param data pointer to the float buffer.
     * \param size number of floats to write.
     * \param os   destination stream.
     */
    inline void writeData (float const* data, std::size_t size, std::ostream& os)
    {
        writeFloatData(data, size, os);
    }

    /**
     * \brief Convenience overload forwarding to writeDoubleData() with the native descriptor.
     * \param data pointer to the double buffer.
     * \param size number of doubles to write.
     * \param os   destination stream.
     */
    inline void writeData (double const* data, std::size_t size, std::ostream& os)
    {
        writeDoubleData(data, size, os);
    }

    /**
     * \brief Convenience overload forwarding to readIntData() with the native descriptor.
     * \param data destination buffer for ints.
     * \param size number of ints to read.
     * \param is   source stream.
     */
    inline void readData (int * data, std::size_t size, std::istream& is)
    {
        readIntData(data, size, is, FPC::NativeIntDescriptor());
    }

    /**
     * \brief Convenience overload forwarding to readLongData() with the native descriptor.
     * \param data destination buffer for longs.
     * \param size number of longs to read.
     * \param is   source stream.
     */
    inline void readData (Long * data, std::size_t size, std::istream& is)
    {
        readLongData(data, size, is, FPC::NativeLongDescriptor());
    }

    /**
     * \brief Convenience overload forwarding to readFloatData() with the native descriptor.
     * \param data destination buffer for floats.
     * \param size number of floats to read.
     * \param is   source stream.
     */
    inline void readData (float * data, std::size_t size, std::istream& is)
    {
        readFloatData(data, size, is, FPC::Native32RealDescriptor());
    }

    /**
     * \brief Convenience overload forwarding to readDoubleData() with the native descriptor.
     * \param data destination buffer for doubles.
     * \param size number of doubles to read.
     * \param is   source stream.
     */
    inline void readData (double * data, std::size_t size, std::istream& is)
    {
        readDoubleData(data, size, is, FPC::Native64RealDescriptor());
    }
}

#endif
