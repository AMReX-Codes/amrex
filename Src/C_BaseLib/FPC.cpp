//BL_COPYRIGHT_NOTICE

//
// $Id: FPC.cpp,v 1.3 2000-04-24 17:52:34 car Exp $
//

#include <FPC.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

//
// FP orders.
//
const int FPC::normal_float_order[]     = { 1, 2, 3, 4 };
const int FPC::reverse_float_order[]    = { 4, 3, 2, 1 };
const int FPC::reverse_float_order_2[]  = { 2, 1, 4, 3 };
const int FPC::normal_double_order[]    = { 1, 2, 3, 4, 5, 6, 7, 8 };
const int FPC::reverse_double_order[]   = { 8, 7, 6, 5, 4, 3, 2, 1 };
const int FPC::reverse_double_order_2[] = { 2, 1, 4, 3, 6, 5, 8, 7 };
const int FPC::cray_float_order[]       = { 1, 2, 3, 4, 5, 6, 7, 8 };

//
// Floating point formats.
//
const long FPC::ieee_float[]  = { 32L,  8L, 23L, 0L, 1L,  9L, 0L,   0x7FL };
const long FPC::ieee_double[] = { 64L, 11L, 52L, 0L, 1L, 12L, 0L,  0x3FFL };
const long FPC::cray_float[]  = { 64L, 15L, 48L, 0L, 1L, 16L, 1L, 0x4000L };


//
// Every copy of the library will have exactly one
// `nativeLongDescriptor' and `nativeRealDescriptor' compiled into it.
// Each machine on which BoxLib runs MUST have them defined below.
//


const
IntDescriptor&
FPC::NativeLongDescriptor ()
{
#if defined(__alpha) || defined(__i486__) || defined(WIN32) || defined(i386) || defined(__i386__)
    static const IntDescriptor nld(sizeof(long), IntDescriptor::ReverseOrder);
#endif

#ifdef _CRAY1
    static const IntDescriptor nld(sizeof(long), IntDescriptor::NormalOrder);
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(_CRAYT3E)  || \
    defined(__hpux)
    static const IntDescriptor  nld(sizeof(long), IntDescriptor::NormalOrder);
#endif

    return nld;
}

const
RealDescriptor&
FPC::NativeRealDescriptor ()
{
#if defined(__alpha) || defined(__i486__) || defined(WIN32) || defined(i386) || defined(__i386__)
#ifdef BL_USE_FLOAT
    static const RealDescriptor nrd(ieee_float, reverse_float_order, 4);
#else
    static const RealDescriptor nrd(ieee_double, reverse_double_order, 8);
#endif
#endif

#ifdef _CRAY1
    static const RealDescriptor nrd(cray_float, cray_float_order, 8);
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(_CRAYT3E)  || \
    defined(__hpux)
#ifdef BL_USE_FLOAT
    static const RealDescriptor nrd(ieee_float, normal_float_order, 4);
#else
    static const RealDescriptor nrd(ieee_double, normal_double_order, 8);
#endif
#endif

    return nrd;
}

const
RealDescriptor&
FPC::CrayRealDescriptor ()
{
    static const RealDescriptor crd(cray_float, cray_float_order, 8);
    return crd;
}

const
RealDescriptor&
FPC::Ieee32NormalRealDescriptor ()
{
    static const RealDescriptor i32rd(ieee_float, normal_float_order, 4);
    return i32rd;
}

const
RealDescriptor&
FPC::Ieee64NormalRealDescriptor ()
{
    static const RealDescriptor i64rd(ieee_double, normal_double_order, 8);
    return i64rd;
}


//
// TODO -- add more machine descriptions.
//
#if !(defined(__alpha)  || \
      defined(_CRAY1)   || \
      defined(_CRAYT3E) || \
      defined(__sgi)    || \
      defined(__sun)    || \
      defined(__i486__) || \
      defined(i386)     || \
      defined(__i386__) || \
      defined(__hpux)   || \
      defined(_MSC_VER) || \
      defined(_AIX))
#error We do not yet support FAB I/O on this machine
#endif

#ifdef BL_NAMESPACE
}
#endif

