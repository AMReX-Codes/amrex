#ifndef BL_IEEE_C_H
#define BL_IEEE_C_H

/***********************************************************************

THIS CODE WAS STOLEN FROM NELSON BEEBE by
CHARLES RENDLEMAN 

This file defines typedefs and symbols for interfacing C IEEE 754
support functions to Fortran code.

The model representation of a t-digit floating-point number is

	x = (-1)**s * 0.d1 d2 d3 ... dt * beta**e

where the digits dk satisfy

	0 <= dk < beta

The fractional part, which we call the significand, is defined to lie
in the range

	1/beta <= significand < 1

For IEEE floating-point, with its hidden bit and denormalized numbers,
we adjust parameters to conform to our model.  Denormalized numbers
are normalized by expanding their exponent range.

IEEE floating point arithmetic has these formats, where s is the sign
bit, e is an exponent bit, and f is a fraction bit.

Single precision:
	seee eeee efff ffff ffff ffff ffff ffff

	significand = 1.fff ffff ffff ffff ffff ffff (1 + 23 bits)

	exponent bias = 127

Double precision:
	seee eeee eeee ffff ffff ffff ffff ffff
	ffff ffff ffff ffff ffff ffff ffff ffff

	significand = 1.ffff ffff ffff ffff ffff ffff ffff ffff
			ffff ffff ffff ffff ffff (1 + 52 bits)

	exponent bias = 1023

Here are some sample IEEE bit patterns:

========================================================================
			    LITTLE ENDIAN

Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	        Infinity	0x7f800000
	       -Infinity	0xff800000
	             NaN	0x7f80ffff

Double precision
	               0	0x00000000 00000000
	               1	0x00000000 3ff00000
	              -1	0x00000000 bff00000
	               2	0x00000000 40000000
	              -2	0x00000000 c0000000
	     1.11022e-16	0x00000002 3ca00000	eps(1.0)
	    -1.11022e-16	0x00000002 bca00000
	    2.22507e-308	0x00000000 00100000	smallest normal
	   -2.22507e-308	0x00000000 80100000
	    4.94066e-324	0x00000001 00000000	smallest subnormal
	   -4.94066e-324	0x00000001 80000000
	    1.79769e+308	0xffffffff 7fefffff	largest normal
	   -1.79769e+308	0xffffffff ffefffff
	        Infinity	0x00000000 7ff00000
	       -Infinity	0x00000000 fff00000
	             NaN	0xffffffff 7ff7ffff

========================================================================
			     BIG ENDIAN
Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	             Inf	0x7f800000
	            -Inf	0xff800000
	             NaN	0x7fffffff

Double precision
	               0	0x00000000 00000000
	               1	0x3ff00000 00000000
	              -1	0xbff00000 00000000
	               2	0x40000000 00000000
	              -2	0xc0000000 00000000
	     1.11022e-16	0x3ca00000 00000002	eps(1.0)
	    -1.11022e-16	0xbca00000 00000002
	    2.22507e-308	0x00100000 00000000	smallest normal
	   -2.22507e-308	0x80100000 00000000
	    4.94066e-324	0x00000000 00000001	smallest subnormal
	   -4.94066e-324	0x80000000 00000001
	    1.79769e+308	0x7fefffff ffffffff	largest normal
	   -1.79769e+308	0xffefffff ffffffff
	             Inf	0x7ff00000 00000000
	            -Inf	0xfff00000 00000000
	             NaN	0x7fffffff ffffffff

========================================================================

***********************************************************************/

/* type mappings from Fortran to C */
typedef double BL_FORT_DOUBLE;
typedef float  BL_FORT_SINGLE;

#if defined(__alpha)
typedef int BL_FORT_INTEGER;	/* need 32-bit integers, not 64-bit ones */
#else
typedef long BL_FORT_INTEGER;
#endif

typedef int BL_FORT_LOGICAL;	/* use int, not long, to avoid conflicts on HP-UX with */
			/* system header file declarations of isinf(), isnan() */

#if defined(BL_FORT_USE_UPPERCASE)
#define sadx	SADX
#define dadx	DADX
#define deps	DEPS
#define seps	SEPS

#define dintxp	DINTXP
#define dsetxp	DSETXP

#define sintxp	SINTXP
#define ssetxp	SSETXP

#define disden	DISDEN
#define disinf	DISINF
#define disnan	DISNAN

#define sisden	SISDEN
#define sisinf	SISINF
#define sisnan	SISNAN

#elif defined(BL_FORT_USE_UNDERSCORE)

#define sadx	sadx_
#define dadx	dadx_
#define deps	deps_
#define seps	seps_

#define dintxp	dintxp_
#define dsetxp	dsetxp_

#define sintxp	sintxp_
#define ssetxp	ssetxp_

#define disden	disden_
#define disinf	disinf_
#define disnan	disnan_

#define sisden	sisden_
#define sisinf	sisinf_
#define sisnan	sisnan_

#elif defined(BL_FORT_USE_DBL_UNDERSCORE)

#define sadx	sadx_
#define dadx	dadx_
#define deps	deps_
#define seps	seps_

#define dintxp	dintxp_
#define dsetxp	dsetxp_

#define sintxp	sintxp_
#define ssetxp	ssetxp_

#define disden	disden_
#define disinf	disinf_
#define disnan	disnan_

#define sisden	sisden_
#define sisinf	sisinf_
#define sisnan	sisnan_

#endif /* (_FTN_NAME_TYPE == 3) */

#if defined(__cplusplus)
extern "C" {
#endif

BL_FORT_DOUBLE dadx(BL_FORT_DOUBLE* x, BL_FORT_INTEGER* n);
BL_FORT_DOUBLE deps(BL_FORT_DOUBLE* x);

BL_FORT_SINGLE sadx(BL_FORT_SINGLE* x, BL_FORT_INTEGER* n);
BL_FORT_SINGLE seps(BL_FORT_SINGLE* x);

BL_FORT_INTEGER dintxp(BL_FORT_DOUBLE* x);
BL_FORT_DOUBLE dsetxp(BL_FORT_DOUBLE* x, BL_FORT_INTEGER* n);

BL_FORT_LOGICAL disden(BL_FORT_DOUBLE* x);
BL_FORT_LOGICAL disinf(BL_FORT_DOUBLE* x);
BL_FORT_LOGICAL disnan(BL_FORT_DOUBLE* x);

BL_FORT_INTEGER sintxp(BL_FORT_SINGLE* x);
BL_FORT_SINGLE ssetxp(BL_FORT_SINGLE* x, BL_FORT_INTEGER* n);

BL_FORT_LOGICAL sisden(BL_FORT_SINGLE* x);
BL_FORT_LOGICAL sisinf(BL_FORT_SINGLE* x);
BL_FORT_LOGICAL sisnan(BL_FORT_SINGLE* x);

#if defined(__cplusplus)
};
#endif

#endif
