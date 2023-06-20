#ifndef AMREX_VISMF_BUFFER_H_
#define AMREX_VISMF_BUFFER_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_Vector.H>

namespace amrex {

class VisMFBuffer
{
public:
    //! We try to do I/O with buffers of this size.
    enum { IO_Buffer_Size = 262144 * 8 };
    //! The type of a char buffer required by [p]setbuf().
#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif
    //! A simple character buffer for setbuf() usage.
    typedef Vector<Setbuf_Char_Type> IO_Buffer;

    static Long GetIOBufferSize () { return ioBufferSize; }
    static void SetIOBufferSize (Long iobuffersize) {
        BL_ASSERT(iobuffersize > 0);
        ioBufferSize = iobuffersize;
    }

protected:

    static AMREX_EXPORT Long ioBufferSize;   //!< ---- the settable buffer size
};

}

#endif
