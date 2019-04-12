
#ifndef AMREX_SPACE_H
#define AMREX_SPACE_H

/* for backward compatibility */
#if defined(BL_SPACEDIM) && !defined(AMREX_SPACEDIM)
#define AMREX_SPACEDIM BL_SPACEDIM
#endif
#if !defined(BL_SPACEDIM)
#define BL_SPACEDIM AMREX_SPACEDIM
#endif

/*
  The AMREX_SPACEDIM macro specifies the dimension of the AMReX library
  to be built.  Acceptable values are either 1, 2, or 3.  This macro
  must be set for every compilation as many of the data structures of
  AMReX rely on it being set.

  Code that depends on the dimensionality of the problem that you're
  trying to solve should be cast in terms of this macro.
*/

#if !defined(AMREX_SPACEDIM)
#error AMREX_SPACEDIM must be defined
#endif

#if  (AMREX_SPACEDIM != 1 && AMREX_SPACEDIM !=2 && AMREX_SPACEDIM != 3)
#error AMREX_SPACEDIM must be either 1, 2, or 3
#endif

/*
  The SpaceDim C++ integer constant specifies the dimension of the
  AMReX library to be built.  It is initialized to the value of
  the macro AMREX_SPACEDIM.  Hence it will have one of the values 1, 2, or 3.

  C++ code that depends on the dimensionality of the problem that
  you're trying to solve should be cast in terms of this constant,
  instead of the AMREX_SPACEDIM macro, if possible, as constants are
  visible to debuggers while macros are not.
*/
#ifdef __cplusplus
namespace amrex
{
    constexpr int SpaceDim = AMREX_SPACEDIM;
}
#endif


#if AMREX_SPACEDIM==1

/*
  The AMREX_D_EXPR(a,b,c) macro expands to a comma expression that will
  evaluate SpaceDim of its arguments, counting from the left; i.e.
  AMREX_SPACEDIM==1 implies only a is evaluated, AMREX_SPACEDIM==2
  implies both a and b, and AMREX_SPACEDIM==3 implies all three
  arguments are evaluated.  This macro can be used to enhance
  the portability of code that works with structures whose size
  depends on AMREX_SPACEDIM.
           
  For example:

     AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);

  More explicitly, this macro evaluates to one of

     ((void)((a),0))

     ((void)((a),(b),0))

     ((void)((a),(b),(c),0))

  corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
*/
#  define AMREX_D_EXPR(a,b,c) ((void)((a),0))

/*
  The AMREX_D_DECL(a,b,c) macro expands to a comma-separated list of
  1, 2, or all 3 of the arguments of the call, depending on the
  value of AMREX_SPACEDIM.  This can be used to write portable
  function calls that depend on 1, 2, or 3 arguments, corresponding
  to the value of AMREX_SPACEDIM.

  For example:

     return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));

  More explicitly, this macro evaluates to one of

     a

     a,b

     a,b,c

  corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
*/
#  define AMREX_D_DECL(a,b,c) a

/*
  The AMREX_D_TERM(a,b,c) macro expands to a whitespace-separated list of
  1, 2, or all 3 of the arguments of the call, depending on the value
  of AMREX_SPACEDIM.  This can be used to write program logic
  that depend on 1, 2, or 3 arguments, corresponding to
  the value of AMREX_SPACEDIM.

  For example:

     return AMREX_D_TERM(len[0], *len[1], *len[2]);

  More explicitly, this macro evaluates to one of

     a

     a b

     a b c

  corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
*/
#  define AMREX_D_TERM(a,b,c) a

/*
  The AMREX_D_PICK(a,b,c) macro expands to a single result equal to the
  1st, 2nd, or 3rd of the arguments of the call, depending on the value
  of AMREX_SPACEDIM.  This can be used to write program logic
  that depend on the value of AMREX_SPACEDIM.

  For example:

     maxsize = AMREX_D_PICK(1024, 128, 32);

  More explicitly, this macro evaluates to one of

     maxsize = 1024;

     maxsize = 128;

     maxsize = 32;

  corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
*/
#  define AMREX_D_PICK(a,b,c) a
#elif AMREX_SPACEDIM==2
#  define AMREX_D_EXPR(a,b,c) ((void)((a),(b),0))
#  define AMREX_D_DECL(a,b,c) a,b
#  define AMREX_D_TERM(a,b,c) a b
#  define AMREX_D_PICK(a,b,c) b
#elif AMREX_SPACEDIM==3
#  define AMREX_D_EXPR(a,b,c) ((void)((a),(b),(c),0))
#  define AMREX_D_DECL(a,b,c) a,b,c
#  define AMREX_D_TERM(a,b,c) a b c
#  define AMREX_D_PICK(a,b,c) c
#endif

#ifndef AMREX_XSDK

/* for backward compatibility */
#if (AMREX_SPACEDIM == 1)
#  define D_DECL(a,b,c) a
#  define D_EXPR(a,b,c) ((void)((a),0))
#  define D_PICK(a,b,c) a
#  define D_TERM(a,b,c) a
#elif (AMREX_SPACEDIM == 2)
#  define D_DECL(a,b,c) a,b
#  define D_EXPR(a,b,c) ((void)((a),(b),0))
#  define D_PICK(a,b,c) b
#  define D_TERM(a,b,c) a b
#elif (AMREX_SPACEDIM == 3)
#  define D_DECL(a,b,c) a,b,c
#  define D_EXPR(a,b,c) ((void)((a),(b),(c),0))
#  define D_PICK(a,b,c) c
#  define D_TERM(a,b,c) a b c
#endif

#endif

#endif /*BL_SPACE_H*/
