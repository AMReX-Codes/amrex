
#ifndef BL_FPC_H
#define BL_FPC_H

#include <AMReX_FabConv.H>

namespace amrex {

/**
* \brief A Collection of Floating-Point Constants Supporting FAB I/O
*
* This class is a poor-person's namespace of floating-point constants
* used in support of FAB I/O.  Since we can't assume the existence of
* namespaces, and we don't like global constants, we make them static
* constant data members of this class.
*/

class FPC
{
public:
    //! 4-element array representing "normal" float order: {1,2,3,4}
    static const int normal_float_order[];

    //! 4-element array representing "reverse" float order: {4,3,2,1}
    static const int reverse_float_order[];

    //! Another "reversed" FP order: {2,1,4,3}
    static const int reverse_float_order_2[];

    //! The "normal" double order: {1,2,3,4,5,6,7,8}
    static const int normal_double_order[];

    //! A "reversed" double order: {8,7,6,5,4,3,2,1}
    static const int reverse_double_order[];

    //! Another "reversed" double order: {2,1,4,3,6,5,8,7}
    static const int reverse_double_order_2[];

    /** \brief Array detailing the format of IEEE 32-bit normal order floats.
    *   In general, here's what the various indices in "format" array means:
    *      format[0] = number of bits per number
    *      format[1] = number of bits in exponent
    *      format[2] = number of bits in mantissa
    *      format[3] = start bit of sign
    *      format[4] = start bit of exponent
    *      format[5] = start bit of mantissa
    *      format[6] = high order mantissa bit (CRAY needs this)
    *      format[7] = bias of exponent
    */
    static const long ieee_float[];

    //! Array detailing the format of IEEE 64-bit normal order doubles.
    static const long ieee_double[];

    /**
    * \brief Returns a constant reference to an IntDescriptor describing
    * the native "long" under which AMReX was compiled.  Each
    * AMReX library will have exactly one of these compiled into it.
    */
    static const IntDescriptor&  NativeLongDescriptor ();

    /**
    * \brief Returns a constant reference to an IntDescriptor describing
    * the native "int" under which AMReX was compiled.  Each
    * AMReX library will have exactly one of these compiled into it.
    */
    static const IntDescriptor&  NativeIntDescriptor ();

    /**
    * \brief Returns a constant reference to a RealDescriptor describing
    * the native Real under which AMReX was compiled.  Each AMReX
    * library will have exactly one of these compiled into it.
    * Note that "Real" will be one of "float" or "double" depending
    * on how the version of AMReX was built.
    */
    static const RealDescriptor& NativeRealDescriptor ();

    /**
    * \brief NativeRealDescriptor is equivalent to Native32RealDescriptor
    * if BL_FLOAT is used. Otherwise, it is equivalent to Native64RealDescriptor.
    */
    static const RealDescriptor& Native32RealDescriptor ();
    static const RealDescriptor& Native64RealDescriptor ();

    /**
    * \brief Returns a constant reference to a RealDescriptor detailing
    * the IEEE 32-bit normal FP format.
    */
    static const RealDescriptor& Ieee32NormalRealDescriptor ();

    /**
    * \brief Returns a constant reference to a RealDescriptor detailing
    * the IEEE 64-bit normal FP format.
    */
    static const RealDescriptor& Ieee64NormalRealDescriptor ();
};

}

#endif /*BL_FPC_H*/
