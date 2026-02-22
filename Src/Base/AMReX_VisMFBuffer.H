#ifndef AMREX_VISMF_BUFFER_H_
#define AMREX_VISMF_BUFFER_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_Vector.H>

/**
 * \file AMReX_VisMFBuffer.H
 * \brief Shared buffer logic for VisMF I/O.
 */

namespace amrex {

class VisMFBuffer
{
public:
    //! We try to do I/O with buffers of this size.
    static constexpr int  IO_Buffer_Size = 262144 * 8;
    //! The type of a char buffer required by [p]setbuf().
#ifdef BL_SETBUF_SIGNED_CHAR
    using Setbuf_Char_Type = signed char;
#else
    using Setbuf_Char_Type = char;
#endif
    //! A simple character buffer for setbuf() usage.
    using IO_Buffer = Vector<Setbuf_Char_Type>;

    //! Return the current global buffer size used by VisMF I/O helpers.
    static Long GetIOBufferSize () { return ioBufferSize; }
    /**
     * \brief Override the default VisMF I/O buffer size (must be positive).
     * \param iobuffersize Desired buffer size in bytes.
     */
    static void SetIOBufferSize (Long iobuffersize) {
        BL_ASSERT(iobuffersize > 0);
        ioBufferSize = iobuffersize;
    }

    /**
     * \brief Release the memory held by a buffer.
     * \param buf Buffer instance to clear and deallocate.
     */
    static void ClearBuffer (IO_Buffer& buf) {
        IO_Buffer().swap(buf);
    }

protected:

    static AMREX_EXPORT Long ioBufferSize;   //!< ---- the settable buffer size
};

}

#endif
