
#include <AMReX_FPC.H>

namespace amrex {

//
// FP orders.
//
const int FPC::normal_float_order[]     = { 1, 2, 3, 4 };
const int FPC::reverse_float_order[]    = { 4, 3, 2, 1 };
const int FPC::reverse_float_order_2[]  = { 2, 1, 4, 3 };
const int FPC::normal_double_order[]    = { 1, 2, 3, 4, 5, 6, 7, 8 };
const int FPC::reverse_double_order[]   = { 8, 7, 6, 5, 4, 3, 2, 1 };
const int FPC::reverse_double_order_2[] = { 2, 1, 4, 3, 6, 5, 8, 7 };
//
// Floating point formats.
//
const long FPC::ieee_float[]  = { 32L,  8L, 23L, 0L, 1L,  9L, 0L,   0x7FL };
const long FPC::ieee_double[] = { 64L, 11L, 52L, 0L, 1L, 12L, 0L,  0x3FFL };
//
// Every copy of the library will have exactly one nativeIntDescriptor,
// nativeLongDescriptor, and nativeRealDescriptor compiled into it.
// Each machine on which AMReX runs MUST have them defined below.
//
const
IntDescriptor&
FPC::NativeIntDescriptor ()
{
#if defined(__i486__) || \
    defined(i386) || \
    defined(__i386__) || \
    defined(__x86_64) || \
    defined(__amd64__) || \
    defined(__LITTLE_ENDIAN__) || \
    defined(__powerpc__) || \
    defined(powerpc)
    static const IntDescriptor nld(sizeof(int), IntDescriptor::ReverseOrder);
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(__ppc__) || \
    defined(__ppc64__) || \
    defined(_SX)   || \
    defined(__hpux)
#if !defined(__LITTLE_ENDIAN__)
    static const IntDescriptor  nld(sizeof(int), IntDescriptor::NormalOrder);
#endif
#endif

    return nld;
}

const
IntDescriptor&
FPC::NativeLongDescriptor ()
{
#if defined(__i486__) || \
    defined(i386) || \
    defined(__i386__) || \
    defined(__x86_64) || \
    defined(__amd64__) || \
    defined(__LITTLE_ENDIAN__) || \
    defined(__powerpc__) || \
    defined(powerpc)
    static const IntDescriptor nld(sizeof(long), IntDescriptor::ReverseOrder);
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(__ppc__) || \
    defined(__ppc64__) || \
    defined(_SX)   || \
    defined(__hpux)
#if !defined(__LITTLE_ENDIAN__)
    static const IntDescriptor  nld(sizeof(long), IntDescriptor::NormalOrder);
#endif
#endif

    return nld;
}

const
RealDescriptor&
FPC::NativeRealDescriptor ()
{
#if defined(__i486__) || \
    defined(i386) || \
    defined(__i386__) || \
    defined(__amd64__) || \
    defined(__LITTLE_ENDIAN__) || \
    defined(__x86_64)
#ifdef BL_USE_FLOAT
    static const RealDescriptor nrd(ieee_float, reverse_float_order, 4);
#else
    static const RealDescriptor nrd(ieee_double, reverse_double_order, 8);
#endif
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(__ppc__) || \
    defined(__ppc64__) || \
    defined(_SX)   || \
    defined(__powerpc__) || \
    defined(powerpc) || \
    defined(__hpux)
#if !defined(__LITTLE_ENDIAN__)
#ifdef BL_USE_FLOAT
    static const RealDescriptor nrd(ieee_float, normal_float_order, 4);
#else
    static const RealDescriptor nrd(ieee_double, normal_double_order, 8);
#endif
#endif
#endif

    return nrd;
}

const
RealDescriptor&
FPC::Native32RealDescriptor ()
{
#if defined(__i486__) || \
    defined(i386) || \
    defined(__i386__) || \
    defined(__amd64__) || \
    defined(__LITTLE_ENDIAN__) || \
    defined(__x86_64)
    static const RealDescriptor n32rd(ieee_float, reverse_float_order, 4);
#endif

#if defined(__sgi) || \
    defined(__sun) || \
    defined(_AIX)  || \
    defined(__ppc__) || \
    defined(__ppc64__) || \
    defined(_SX)   || \
    defined(__powerpc__) || \
    defined(powerpc) || \
    defined(__hpux)
#if !defined(__LITTLE_ENDIAN__)
    static const RealDescriptor n32rd(ieee_float, normal_float_order, 4);
#endif
#endif
    return n32rd;
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


#if !(defined(_SX)      || \
      defined(__sgi)    || \
      defined(__sun)    || \
      defined(__i486__) || \
      defined(i386)     || \
      defined(__ppc__) || \
      defined(__ppc64__) || \
      defined(__i386__) || \
      defined(__amd64__) || \
      defined(__x86_64) || \
      defined(__hpux)   || \
      defined(__powerpc__) || \
      defined(powerpc)  || \
      defined(__LITTLE_ENDIAN__)  || \
      defined(_MSC_VER) || \
      defined(_AIX))
#error We do not yet support FAB I/O on this machine
#endif

}
