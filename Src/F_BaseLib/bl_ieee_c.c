#include "bl_ieee_c.h"

/* Fortran BL_FORT_LOGICAL values -- compiler dependent!  Most (all?) UNIX */
/* Fortran compilers use 0 for .FALSE. and non-zero for .TRUE., like C. */

#ifdef _FALSE_
#undef _FALSE_
#endif

#ifdef _TRUE_
#undef _TRUE_
#endif

#define _FALSE_				((BL_FORT_LOGICAL)0)
#define _TRUE_				((BL_FORT_LOGICAL)1)

#define BASE				2 /* number base */

/* stored significand bits in single and double precision */

#define	T_SP				23
#define	T_DP				52

#define BASE_TO_THE_T_SP		((BL_FORT_SINGLE)8388608.0)
#define BASE_TO_THE_T_DP		((BL_FORT_DOUBLE)4503599627370496.0)

#define _SHIFTED_EXPONENT_MASK_SP	0xff
#define EXPONENT_MASK_SP		0x7f800000L
#define EXPONENT_MASK_DP		0x7ff00000L

#define _SHIFTED_EXPONENT_MASK_DP	0x7ff
#define _SIGNIFICAND_MASK_SP		0x007fffffL
#define _SIGNIFICAND_MASK_DP		0x000fffffL

/* Exponent biases such that significand lies in (1/beta) <= significand < 1.
These are 1 less than the IEEE biases, because its stored significand
lies in 1 <= significand < beta due to the hidden bit.  We define them
with a leading underscore because they are for internal use only. */

#define _BIAS_SP			126
#define _BIAS_DP			1022

/* Indexes into two-word BL_FORT_INTEGER array to account for addressing order */

#if i386 || sun386 || __i386__ || __sun386__ || msdos || MIPSEL || __alpha || WIN32
					/* Intel 80xxx or MIPS little endian */
#define DP_LOW				0
#define DP_HIGH				1
#else				/* big endian (MIPS, Motorola, SPARC, ...) */
#define DP_LOW				1
#define DP_HIGH				0
#endif

/* macros to extract (high-order) significand and exponent as integer values */

#define GET_EXPONENT_SP(x) ((((x) >> T_SP) & \
				_SHIFTED_EXPONENT_MASK_SP) - _BIAS_SP)
#define GET_EXPONENT_DP(x) ((((x) >> (T_DP - 32)) & \
				_SHIFTED_EXPONENT_MASK_DP) - _BIAS_DP)

#define SET_EXPONENT_SP(x)	(((x) + _BIAS_SP) << T_SP)
#define SET_EXPONENT_DP(x)	(((x) + _BIAS_DP) << (T_DP - 32))

#define EXPONENT_DENORM_SP		(-_BIAS_SP)
#define EXPONENT_DENORM_DP		(-_BIAS_DP)

#define EXPONENT_INFNAN_SP		(255 - _BIAS_SP)
#define EXPONENT_INFNAN_DP		(2047 - _BIAS_DP)

#define SIGNIFICAND_SP(x)		(((x) & _SIGNIFICAND_MASK_SP))
#define SIGNIFICAND_DP(x)		(((x) & _SIGNIFICAND_MASK_DP))

#define MAX_NORMAL_SP			0x7f7fffffL
#define MAX_NORMAL_DP			0x7fefffffL
#define MAX_NORMAL_Low_DP		0xffffffffL

#define MIN_NORMAL_SP			0x00800000L
#define MIN_DENORMAL_SP			0x00000001L

#define MIN_NORMAL_DP			0x00100000L
#define MIN_NORMAL_Low_DP		0x00000000L
#define MIN_DENORMAL_DP			0x00000000L
#define MIN_DENORMAL_Low_DP		0x00000001L

#define Inf_SP				0x7f800000L
#define NegInf_SP			0xff800000L
#define NaN_SP				0x7fffffffL /* significand is */
						    /* arbitrary non-zero */

/* High-order words for double-precision Infinity and NaN. */
#define Inf_DP				0x7ff00000L
#define Inf_Low_DP			0x00000000L
#define NegInf_DP			0xfff00000L
#define NegInf_Low_DP			0x00000000L
#define NaN_DP				0x7fffffffL /* significand is */
#define NaN_Low_DP			0xffffffffL /* arbitrary non-zero */

#define ISNEG_SP(x)			((x) & 0x80000000L)
#define ISNEG_DP(x)			((x) & 0x80000000L)

#if STDC
#include <stdlib.h>
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif


BL_FORT_DOUBLE
dadx(BL_FORT_DOUBLE *x, BL_FORT_INTEGER *n)
{
    BL_FORT_INTEGER delta_e;
    BL_FORT_INTEGER e;
    BL_FORT_DOUBLE_PARTS w;
    BL_FORT_DOUBLE new_x;

    w.r = *x;
    e = GET_EXPONENT_DP(w.i[DP_HIGH]);		/* extract old exponent */

    if (*x == 0.0)			/* zero remains unchanged */
	/* NO-OP */;
    else if (e == EXPONENT_INFNAN_DP)	/* original number is NaN or Inf */
	/* NO-OP */;			/* so no change in value */
    else if (e == EXPONENT_DENORM_DP)	/* original number is denormalized */
    {
	delta_e = *n - T_DP;
	new_x = *x * BASE_TO_THE_T_DP;	/* scale into normal range */
	w.r = dadx(&new_x, &delta_e);
    }
    else				/* original number is normal */
    {
	e += *n;			/* new exponent */
	w.r = dsetxp(x,&e);
    }
    return (w.r);
}


BL_FORT_DOUBLE
deps(BL_FORT_DOUBLE *x)
{
    BL_FORT_DOUBLE		half = 0.5;
    BL_FORT_INTEGER			n = dintxp(x) - T_DP;
    BL_FORT_INTEGER			zero = 0;
    BL_FORT_DOUBLE_PARTS	w;

    if ( (*x == 0.0) || disden(x) )
    {
	w.i[DP_HIGH] = MIN_DENORMAL_DP;
	w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	return (w.r);
    }
    else if ( disnan(x) || disinf(x) )
	return *x;
    else if ( (*x < 0.0) && (dsetxp(x,&zero) == -half) )
    {					/* special boundary case */
	w.r = dsetxp(&half,(n--, &n));
	if (w.r == 0.0)		/* then x = 2.22507e-308 0x00000000 00100000 */
	{
	    w.i[DP_HIGH] = MIN_DENORMAL_DP; /* and setxp() -> 0, so fix up */
	    w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	}
	return (w.r);
    }
    else
	return (dsetxp(&half,&n));
}

BL_FORT_LOGICAL
disinf(BL_FORT_DOUBLE *x)
{
    BL_FORT_DOUBLE_PARTS w;

    w.r = *x;
    return (((GET_EXPONENT_DP(w.i[DP_HIGH]) == EXPONENT_INFNAN_DP) &&
	     (SIGNIFICAND_DP(w.i[DP_HIGH]) == 0) && (w.i[DP_LOW] == 0))
	    ? _TRUE_ : _FALSE_);
}

BL_FORT_LOGICAL
disden(BL_FORT_DOUBLE *x)
{
    return (dintxp(x) <= EXPONENT_DENORM_DP) ? _TRUE_ : _FALSE_;
}

BL_FORT_LOGICAL
disnan(BL_FORT_DOUBLE *x)
{
#if defined(__alpha)
    if (disden(x))
	return (_FALSE_);
    else
	return ((*x != *x) ? _TRUE_ : _FALSE_);
#else
    return ((*x != *x) ? _TRUE_ : _FALSE_);
#endif
}

BL_FORT_DOUBLE
dsetxp(BL_FORT_DOUBLE *x, BL_FORT_INTEGER *n)
{
    BL_FORT_DOUBLE_PARTS w;

    if (disden(x))			/* then x is denormalized */
        w.r = *x * BASE_TO_THE_T_DP;	/* so must first normalize fraction */
    else
	w.r = *x;			/* x is normal, NaN, or Inf */
    if (w.r == 0.0)			/* zeroes must be left intact */
        /* NO-OP */;
    else if (disnan(x))		/* NaNs must be left intact */
        /* NO-OP */;
    else if (disinf(x))		/* Infs must be left intact */
        /* NO-OP */;
    else if (*n <= (EXPONENT_DENORM_DP - T_DP)) /* underflow to zero */
	w.r = (BL_FORT_DOUBLE)0.0;
    else if (*n <= EXPONENT_DENORM_DP)	/* generate denormalized value */
    {
	w.i[DP_HIGH] &= ~EXPONENT_MASK_DP; /* clear exponent field */
	w.i[DP_HIGH] |= SET_EXPONENT_DP(*n + T_DP); /* set new */
					/* exponent of scaled normal value */
	w.r /= BASE_TO_THE_T_DP;	/* and denormalize by division */
    }
    else if (*n >= EXPONENT_INFNAN_DP)	/* generate infinity */
    {
	w.i[DP_HIGH] = ISNEG_DP(w.i[DP_HIGH]) ? NegInf_DP : Inf_DP;
	w.i[DP_LOW] = ISNEG_DP(w.i[DP_HIGH]) ? NegInf_Low_DP : Inf_Low_DP;
    }
    else				/* result is normal number */
    {					/* and exponent in range */
	w.i[DP_HIGH] &= ~EXPONENT_MASK_DP;	/* clear exponent field */
	w.i[DP_HIGH] |= SET_EXPONENT_DP(*n);	/* set new exponent */
    }

    return (w.r);
}

BL_FORT_INTEGER
dintxp(BL_FORT_DOUBLE *x)
{
    BL_FORT_DOUBLE_PARTS w;
    register BL_FORT_INTEGER e;

    w.r = *x;
    e = GET_EXPONENT_DP(w.i[DP_HIGH]);

    if (*x == 0.0)			/* handle zero specially */
	e = 0;
    else if (e == EXPONENT_DENORM_DP)	/* have denormalized number */
    {
	w.r *= BASE_TO_THE_T_DP;	/* make normal number */
	e = GET_EXPONENT_DP(w.i[DP_HIGH]) - T_DP;
    }
    return (e);
}

BL_FORT_SINGLE
sadx(BL_FORT_SINGLE *x, BL_FORT_INTEGER *n)
{
    BL_FORT_INTEGER delta_e;
    BL_FORT_INTEGER e;
    BL_FORT_SINGLE_PARTS w;
    BL_FORT_SINGLE new_x;

    w.r = *x;
    e = GET_EXPONENT_SP(w.i);		/* extract old exponent */

    if (*x == 0.0)			/* zero remains unchanged */
	/* NO-OP */;
    else if (e == EXPONENT_INFNAN_SP)	/* original number is NaN or Inf */
	/* NO-OP */;			/* so no change in value */
    else if (e == EXPONENT_DENORM_SP)	/* original number is denormalized */
    {
	delta_e = *n - T_SP;
	new_x = *x * BASE_TO_THE_T_SP;	/* scale into normal range */
	w.r = sadx(&new_x, &delta_e);
    }
    else				/* original number is normal */
    {
	e += *n;			/* new exponent */
	w.r = ssetxp(x,&e);
    }
    return (w.r);
}

BL_FORT_SINGLE
seps(BL_FORT_SINGLE *x)
{
    BL_FORT_SINGLE	half = 0.5;
    BL_FORT_INTEGER	n = sintxp(x) - T_SP;
    BL_FORT_INTEGER	zero = 0;
    BL_FORT_SINGLE_PARTS	w;

    if ( (*x == 0.0) || sisden(x) )
    {
	w.i = MIN_DENORMAL_SP;
	return (w.r);
    }
    else if ( sisnan(x) || sisinf(x) )
	return *x;
    else if ( (*x < 0.0) && (ssetxp(x,&zero) == -half) )
    {					/* special boundary case */
	w.r = ssetxp(&half,(n--, &n));
	if (w.r == 0.0)			/* then x = -1.17549e-38 80800000 */
	    w.i = MIN_DENORMAL_SP;	/* and setxp() -> 0, so fix up */
	return (w.r);
    }
    else
	return (ssetxp(&half,&n));
}

BL_FORT_LOGICAL
sisinf(BL_FORT_SINGLE *x)
{
    BL_FORT_SINGLE_PARTS w;

    w.r = *x;
    return (((GET_EXPONENT_SP(w.i) == EXPONENT_INFNAN_SP) &&
	     (SIGNIFICAND_SP(w.i) == 0))
	    ? _TRUE_ : _FALSE_);
}


BL_FORT_LOGICAL
sisden(BL_FORT_SINGLE *x)
{
    return (sintxp(x) <= EXPONENT_DENORM_SP) ? _TRUE_ : _FALSE_;
}

BL_FORT_LOGICAL
sisnan(BL_FORT_SINGLE *x)
{
#if defined(__alpha)
    if ( sisden(x) )
	return (_FALSE_);
    else
	return ((*x != *x) ? _TRUE_ : _FALSE_);
#else
    return ((*x != *x) ? _TRUE_ : _FALSE_);
#endif
}

BL_FORT_SINGLE
ssetxp(BL_FORT_SINGLE *x, BL_FORT_INTEGER *n)
{
    BL_FORT_SINGLE_PARTS w;

    if ( sisden(x) )			/* then x is denormalized */
        w.r = *x * BASE_TO_THE_T_SP;	/* so must first normalize fraction */
    else
	w.r = *x;			/* x is normal, NaN, or Inf */
    if (w.r == 0.0)			/* zeroes must be left intact */
        /* NO-OP */;
    else if (sisnan(x))			/* NaNs must be left intact */
        /* NO-OP */;
    else if (sisinf(x))			/* Infs must be left intact */
        /* NO-OP */;
    else if (*n <= (EXPONENT_DENORM_SP - T_SP)) /* underflow to zero */
	w.r = (BL_FORT_SINGLE)0.0;
    else if (*n <= EXPONENT_DENORM_SP)	/* generate denormalized value */
    {
	w.i &= ~EXPONENT_MASK_SP;	/* clear exponent field */
	w.i |= SET_EXPONENT_SP(*n + T_SP);	/* set new exponent */
					/* of scaled normal value */
	w.r /= BASE_TO_THE_T_SP;	/* and denormalize by division */
    }
    else if (*n >= EXPONENT_INFNAN_SP)	/* generate infinity */
	w.i = ISNEG_SP(w.i) ? NegInf_SP : Inf_SP;
    else				/* result is normal number */
    {					/* and exponent in range */
	w.i &= ~EXPONENT_MASK_SP;	/* clear exponent field */
	w.i |= SET_EXPONENT_SP(*n);	/* set new exponent */
    }

    return (w.r);
}


BL_FORT_INTEGER
sintxp(BL_FORT_SINGLE *x)
{
    BL_FORT_SINGLE_PARTS w;
    register BL_FORT_INTEGER e;

    w.r = *x;
    e = GET_EXPONENT_SP(w.i);

    if (*x == 0.0)			/* handle zero specially */
	e = 0;
    else if (e == EXPONENT_DENORM_SP)	/* have denormalized number */
    {
	w.r *= BASE_TO_THE_T_SP;	/* make normal number */
	e = GET_EXPONENT_SP(w.i) - T_SP;
    }
    return (e);
}
