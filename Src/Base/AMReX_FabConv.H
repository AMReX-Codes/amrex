
#ifndef BL_FABCONV_H
#define BL_FABCONV_H

#include <iosfwd>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>

namespace amrex {

//
// A Descriptor of the Long Integer type

/**
* This class is meant to hold all information needed to completely
* describe the "int" or "long" type on a machine.  To describe an integer both
* the number of bytes and their ordering, relative to canonical
* ordering 1 .. sizeof(long), needs to be specified.
* This allows us to write out integers in the native format on a machine,
* and then by also saving the IntDescriptor, we can read them back in on
* another machine and have enough information to construct the exact same
* values.
*/

class IntDescriptor
{

public:
    /**
    * \brief An enumeration describing the two orderings
    * that we currently support: NormalOrder and ReverseOrder.
    * Other orderings may be added as AMReX is ported to run on
    * machines with non-standard orderings.
    */
    enum Ordering { NormalOrder = 1, ReverseOrder = 2 };
    /**
    * \brief The default constructor.  Does not build a proper IntDescriptor.
    * This should only be used when you need to build a generic
    * IntDescriptor in order to be able to read in a specific
    * one from an istream.
    */
    IntDescriptor ();
    //! Construct a specific IntDescriptor.
    IntDescriptor (long     nb,
                   Ordering ordering = NormalOrder);
    //! Returns the ordering.
    Ordering order () const;
    //! Returns the number of bytes.
    int numBytes () const;
    //! The equality operator.
    bool operator== (const IntDescriptor& id) const;
    //! The inequality operator.
    bool operator!= (const IntDescriptor& id) const;

private:
    friend std::istream& operator>> (std::istream& is, IntDescriptor& id);
    long     numbytes;
    Ordering ord;
};

//!
//! Write out an IntDescriptor to an ostream in ASCII.
//!
std::ostream& operator<< (std::ostream& os, const IntDescriptor& id);
//!
//! Read in an IntDescriptor from an istream.
//!
std::istream& operator>> (std::istream& is, IntDescriptor& id);


  //A Descriptor of the Real Type

/**
* \brief  This class is meant to hold all information needed to completely
*  describe the "Real" floating-point type on a machine.  By "Real" here we
*  mean either the "float" or "double" type that this version of AMReX
*  was built with, which corresponds to whether BL_USE_FLOAT or
*  BL_USE_DOUBLE was used to build the version of the library.
*
*  To describe a "Real" type two arrays are needed: one detailing the ordering
*  of the bytes in the Real, relative to the canonical ordering
*  1 .. sizeof(Real) and the other detailing the format of the floating-point
*  number.
*
*  The array detailing the format of a floating-point number is an eight-element
*  array of longs containing the following information:
*
*          format[0] = number of bits per number
*          format[1] = number of bits in exponent
*          format[2] = number of bits in mantissa
*          format[3] = start bit of sign
*          format[4] = start bit of exponent
*          format[5] = start bit of mantissa
*          format[6] = high order mantissa bit (CRAY needs this)
*          format[7] = bias of exponent
*
*  This allows us to write out "Real"s in the native format on a machine,
*  and then by also saving the IntDescriptor, we can read them back in on
*  another machine and have enough information to construct the exact same
*  "Real" values, provided the Reals have the same size on the two machines.
*/

class RealDescriptor
{
public:
    /**
    * \brief The default constructor.  Does not build a proper
    * RealDescriptor. This should only be used when you need to
    * build a generic RealDescriptor in order to be able to read
    * in a specific one from an istream.
    */
    RealDescriptor ();
    /**
    * \brief Construct a specific RealDescriptor, passing in the format
    * of the Real, the order of the Real, and the length of the
    * array detailing the ordering.
    */
    RealDescriptor (const long* format,
                    const int*  order,
                    int         order_length);
    //! Returns the format array as a const long*
    const long* format () const&;

    const long* format () && = delete;
    //! Returns const Vector<long> reference to the format array.

    const Vector<long>& formatarray () const&;
    const Vector<long>& formatarray () && = delete;

    //! Returns the order array as a const int*
    const int* order () const&;
    const int* order () && = delete;

    //! Returns const Vector<int> reference to the order array.
    const Vector<int>& orderarray () const&;

    const Vector<int>& orderarray () && = delete;
    //! Returns the number of bytes in the Real.
    int numBytes () const;

    //! The equality operator.
    bool operator== (const RealDescriptor& rd) const;

    //! The inequality operator.
    bool operator != (const RealDescriptor& rd) const;

    //! Set to always fix denormals when converting to native format.
    static void SetFixDenormals ();

    //! Set read and write buffer sizes
    static void SetReadBufferSize (int rbs);
    static void SetWriteBufferSize (int wbs);

    /**
    * \brief Returns a copy of this RealDescriptor on the heap.
    * The user is responsible for deletion.
    */
    RealDescriptor* clone () const;

    /**
    * \brief Returns a RealDescriptor on the heap detailing requested
    * floating-point type.  Here format, precision, and ordering
    * correspond to the enumerations in FABio.  This is here to
    * support reading "old" FABs.  Do NOT use it in new code.
    */
    static RealDescriptor* newRealDescriptor (int         format,
                                              int         precision,
                                              const char* systype,
                                              int         ordering);
    /**
    * \brief Convert nitems in RealDescriptor format to native Real format.
    * The out array is assumed to be large enough to hold the
    * resulting output.
    */
    static void convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       void*                 in,
                                       const RealDescriptor& id);

    /**
    * \brief Read nitems from istream in RealDescriptor format and
    * convert them to the native Real format.  The out array is
    * assumed to be large enough to hold the resulting output.
    */
    static void convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       std::istream&         is,
                                       const RealDescriptor& id);

    /**
    * \brief Convert nitems Reals in native format to RealDescriptor format
    * and write them to the ostream.
    */
    static void convertFromNativeFormat (std::ostream&         os,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od);
    /**
    * \brief Convert nitems Reals in native format to RealDescriptor format.
    * The out array is assumed to be large enough to hold the
    * resulting output.
    */
    static void convertFromNativeFormat (void*                 out,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od);

    /**
    * \brief Convert nitems floats in native format to RealDescriptor format
    * and write them to the ostream.
    */
    static void convertFromNativeFloatFormat (std::ostream&         os,
                                              long                  nitems,
                                              const float*          in,
                                              const RealDescriptor& od);

    /**
    * \brief Convert nitems doubles in native format to RealDescriptor format
    * and write them to the ostream.
    */
    static void convertFromNativeDoubleFormat (std::ostream&         os,
                                               long                  nitems,
                                               const double*         in,
                                               const RealDescriptor& od);

    /**
    * \brief Read nitems from istream in RealDescriptor format and
    * convert them to the native float format.  The out array is
    * assumed to be large enough to hold the resulting output.
    */
    static void convertToNativeFloatFormat (float*                out,
                                            long                  nitems,
                                            std::istream&         is,
                                            const RealDescriptor& id);

    /**
    * \brief Read nitems from istream in RealDescriptor format and
    * convert them to the native double format.  The out array is
    * assumed to be large enough to hold the resulting output.
    */
    static void convertToNativeDoubleFormat (double*               out,
                                             long                  nitems,
                                             std::istream&         is,
                                             const RealDescriptor& id);

private:

    Vector<long> fr;
    Vector<int>  ord;
    static bool bAlwaysFixDenormals;
    static int writeBufferSize;
    static int readBufferSize;
};

//!
//! Write out an RealDescriptor to an ostream in ASCII.
//!
std::ostream& operator<< (std::ostream& os, const amrex::RealDescriptor& id);
//!
//! Read in a RealDescriptor from an istream.
//!
std::istream& operator>> (std::istream& is, amrex::RealDescriptor& id);

}

#endif /*BL_FABCONV_H*/
