//BL_COPYRIGHT_NOTICE

#ifndef BL_FABCONV_H
#define BL_FABCONV_H

//
// $Id: FabConv.H,v 1.2 1997-09-17 22:51:16 lijewski Exp $
//
#include <iostream.h>

#include <Array.H>
#include <Assert.H>
#include <REAL.H>

//
//@Man:
//@Memo: A Descriptor of the Long Integer type
/*@Doc:

  This class is meant to hold all information needed to completely
  describe the "long" type on a machine.  To describe a "long" both
  the number of bytes in the long and their ordering, relative to canonical
  ordering 1 .. sizeof(long), needs to be specified.

  This allows us to write out "long"s in the native format on a machine,
  and then by also saving the IntDescriptor, we can read them back in on
  another machine and have enough information to construct the exact same
  "long" values.
*/

class IntDescriptor
{

public:
    /*@ManDoc: An enumeration describing the two ordering of a "long"
               that we currently support: NormalOrder and ReverseOrder.
               Other ordering may be added as BoxLib is ported to run on
               machines with non-standard "long" orderings.
    */
    enum Ordering { NormalOrder = 1, ReverseOrder = 2 };

    /*@ManDoc: The default constructor.  Does not build a proper IntDescriptor.
               This should only be used when you need to build a generic
               IntDescriptor in order to be able to read in a specific
               one from an istream.
    */
    IntDescriptor ();
    //
    //@ManDoc: Construct a specific IntDescriptor.
    //
    IntDescriptor (long     nb,
                   Ordering ordering = NormalOrder);
    //
    //@ManDoc: Returns the ordering of the "long".
    //
    Ordering order () const;
    //
    //@ManDoc: Returns the number of bytes in a "long".
    //
    int numBytes () const;
    //
    //@ManDoc: The equality operator.
    //
    bool operator== (const IntDescriptor& id) const;
    //
    //@ManDoc: The inequality operator.
    //
    bool operator!= (const IntDescriptor& id) const;
    //
    //@ManDoc: Write out an IntDescriptor to an ostream in ASCII.
    //
    friend ostream& operator<< (ostream&             os,
                                const IntDescriptor& id);
    //
    //@ManDoc: Read in an IntDescriptor from an istream.
    //
    friend istream& operator>> (istream&       is,
                                IntDescriptor& id);
protected:
    long     numbytes;
    Ordering ord;
};

//
//@Man:
//@Memo: A Descriptor of the Real Type
/*@Doc:

  This class is meant to hold all information needed to completely
  describe the "Real" floating-point type on a machine.  By "Real" here we
  mean either the "float" or "double" type that this version of BoxLib
  was built with, which corresponds to whether BL\_USE\_FLOAT or
  BL\_USE\_DOUBLE was used to build the version of the library.
  
  To describe a "Real" type two arrays are needed: one detailing the ordering
  of the bytes in the Real, relative to the canonical ordering
  1 .. sizeof(Real) and the other detailing the format of the floating-point
  number.

  The array detailing the format of a floating-point number is an eight-element
  array of longs containing the following information:

          format[0] = number of bits per number
          format[1] = number of bits in exponent
          format[2] = number of bits in mantissa
          format[3] = start bit of sign
          format[4] = start bit of exponent
          format[5] = start bit of mantissa
          format[6] = high order mantissa bit (CRAY needs this)
          format[7] = bias of exponent

  This allows us to write out "Real"s in the native format on a machine,
  and then by also saving the IntDescriptor, we can read them back in on
  another machine and have enough information to construct the exact same
  "Real" values, provided the Reals have the same size on the two machines.
*/

class RealDescriptor
{
public:
    /*@ManDoc: The default constructor.  Does not build a proper
               RealDescriptor. This should only be used when you need to
               build a generic RealDescriptor in order to be able to read
               in a specific one from an istream.
    */
    RealDescriptor ();

    /*@ManDoc: Construct a specific RealDescriptor, passing in the format
               of the Real, the order of the Real, and the length of the
               array detailing the ordering.
    */
    RealDescriptor (const long* format,
                    const int*  order,
                    int         order_length);
    //
    //@ManDoc: The copy constructor.
    //
    RealDescriptor (const RealDescriptor& rhs);
    //
    //@ManDoc: The assignment operator.
    //
    RealDescriptor& operator= (const RealDescriptor& rhs);
    //
    //@ManDoc: Returns the format array as a `const long*'
    //
    const long* format () const;
    //
    //@ManDoc: Returns const Array<long> reference to the format array.
    //
    const Array<long>& formatarray () const;
    //
    //@ManDoc: Returns the order array as a `const int*'
    //
    const int* order () const;
    //
    //@ManDoc: Returns const Array<int> reference to the order array.
    //
    const Array<int>& orderarray () const;
    //
    //@ManDoc: Returns the number of bytes in the Real.
    //
    int numBytes () const;
    //
    //@ManDoc: The equality operator.
    //
    bool operator== (const RealDescriptor& rd) const;
    //
    //@ManDoc: The inequality operator.
    //
    bool operator != (const RealDescriptor& rd) const;

    /*@ManDoc: Returns a copy of this RealDescriptor on the heap.
               The user is responsible for deletion.
    */
    virtual RealDescriptor* clone () const;

    /*@ManDoc: Returns a RealDescriptor on the heap detailing requested
               floating-point type.  Here format, precision, and ordering
               correspond to the enumerations in FABio.  This is here to
               support reading "old" FABs.  Do NOT use it in new code.
    */
    static RealDescriptor* newRealDescriptor (int         format,
                                              int         precision,
                                              const char* systype,
                                              int         ordering);

    /*@ManDoc: Convert nitems in RealDescriptor format to native Real format.
               The out array is assumed to be large enough to hold the
               resulting output.
    */
    static void convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       void*                 in,
                                       const RealDescriptor& id);

    /*@ManDoc: Read nitems from istream in ReadDescriptor format and
               convert them to the native Real format.  The out array is
               assumed to be large enough to hold the resulting output.
    */
    static void convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       istream&              is,
                                       const RealDescriptor& id);

    /*@ManDoc: Convert nitems Reals in native format to RealDescriptor format
               and write them to the ostream.
    */
    static void convertFromNativeFormat (ostream&              os,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od);

    /*@ManDoc: Convert nitems Reals in native format to RealDescriptor format.
               The out array is assumed to be large enough to hold the
               resulting output.
    */
    static void convertFromNativeFormat (void*                 out,
                                         long                  nitems,
                                         Real*                 in,
                                         const RealDescriptor& od);
    //
    //@ManDoc: Write out an RealDescriptor to an ostream in ASCII.
    //
    friend ostream& operator<< (ostream&              os,
                                const RealDescriptor& id);
    //
    //@ManDoc: Read in a RealDescriptor from an istream.
    //
    friend istream& operator>> (istream&        is,
                                RealDescriptor& id);
protected:
    Array<long> fr;
    Array<int>  ord;
};

//
// Inlines.
//

inline
IntDescriptor::IntDescriptor ()
{}

inline
IntDescriptor::IntDescriptor (long     nb,
                              Ordering o)
    : numbytes(nb),
      ord(o)
{}

inline
IntDescriptor::Ordering
IntDescriptor::order () const
{
    return ord;
}

inline
int
IntDescriptor::numBytes () const
{
    return numbytes;
}

inline
bool
IntDescriptor::operator== (const IntDescriptor& id) const
{
    return ord == id.ord && numbytes == id.numbytes;
}

inline
bool
IntDescriptor::operator!= (const IntDescriptor& id) const
{
    return !operator==(id);
}

inline
RealDescriptor::RealDescriptor ()
{}

inline
RealDescriptor::RealDescriptor (const long* fr_,
                                const int*  ord_,
                                int         ordl_)
    : fr(fr_, 8),
      ord(ord_, ordl_)
{}

inline
RealDescriptor::RealDescriptor (const RealDescriptor& rhs)
    : fr(rhs.fr),
      ord(rhs.ord)
{}

inline
RealDescriptor&
RealDescriptor::operator= (const RealDescriptor& rhs)
{
    fr  = rhs.fr;
    ord = rhs.ord;
    return *this;
}

inline
const long*
RealDescriptor::format () const
{
    assert(fr.length() != 0);
    return fr.dataPtr();
}

inline
const Array<long>&
RealDescriptor::formatarray () const
{
    assert(fr.length() != 0);
    return fr;
}

inline
const int*
RealDescriptor::order () const
{
    assert(ord.length() != 0);
    return ord.dataPtr();
}

inline
const Array<int>&
RealDescriptor::orderarray () const
{
    assert(ord.length() != 0);
    return ord;
}

inline
int
RealDescriptor::numBytes () const
{
    assert(fr.length() != 0);
    return (fr[0] + 7 ) >> 3;
}

inline
bool
RealDescriptor::operator== (const RealDescriptor& rd) const
{
    return fr == rd.fr && ord == rd.ord;
}

inline
bool
RealDescriptor::operator != (const RealDescriptor& rd) const
{
    return !operator==(rd);
}

#endif /*BL_FABCONV_H*/
