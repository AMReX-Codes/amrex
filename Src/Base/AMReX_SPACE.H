/**
 * \file AMReX_SPACE.H
 * \brief Compile-time dimensionality helpers used throughout AMReX.
 */

#ifndef AMREX_SPACE_H
#define AMREX_SPACE_H
#include <AMReX_Config.H>

#if !defined(AMREX_SPACEDIM) && (defined(DOXYGEN) || defined(__DOXYGEN__))
#define AMREX_SPACEDIM 3
#endif

/* for backward compatibility */
#if defined(BL_SPACEDIM) && !defined(AMREX_SPACEDIM)
#define AMREX_SPACEDIM BL_SPACEDIM
#endif
#if !defined(BL_SPACEDIM)
#define BL_SPACEDIM AMREX_SPACEDIM
#endif

/**
 * \brief Dimension macro (`AMREX_SPACEDIM`) that selects 1D, 2D, or 3D builds.
 *
 * The macro must be defined (usually by the build system) for every compilation,
 * since many AMReX data structures depend on it. Code that is dimension-aware
 * should refer to `AMREX_SPACEDIM` (or \ref amrex::SpaceDim) instead of hard-coding values.
 */

#if !defined(AMREX_SPACEDIM)
#error AMREX_SPACEDIM must be defined
#endif

#if  (AMREX_SPACEDIM != 1 && AMREX_SPACEDIM !=2 && AMREX_SPACEDIM != 3)
#error AMREX_SPACEDIM must be either 1, 2, or 3
#endif

/**
 * \brief Dimension constant available in the \c amrex namespace.
 *
 * The SpaceDim C++ integer constant specifies the dimension of the
 * AMReX library to be built.  It is initialized to the value of
 * the macro AMREX_SPACEDIM.  Hence it will have one of the values 1, 2, or 3.
 *
 * C++ code that depends on the dimensionality of the problem that
 * you're trying to solve should be cast in terms of this constant,
 * instead of the AMREX_SPACEDIM macro, if possible, as constants are
 * visible to debuggers while macros are not.
 */
#ifdef __cplusplus
namespace amrex
{
    constexpr int SpaceDim = AMREX_SPACEDIM;
}
#endif


#if AMREX_SPACEDIM==1

/**
 * The AMREX_D_EXPR(a,b,c) macro expands to a comma expression that will
 * evaluate SpaceDim of its arguments, counting from the left; i.e.
 * AMREX_SPACEDIM==1 implies only a is evaluated, AMREX_SPACEDIM==2
 * implies both a and b, and AMREX_SPACEDIM==3 implies all three
 * arguments are evaluated.  This macro can be used to enhance
 * the portability of code that works with structures whose size
 * depends on AMREX_SPACEDIM.
 *
 * For example:
 * \code{.cpp}
 *    AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
 * \endcode
 *
 * More explicitly, this macro evaluates to one of
 *
 * \code{.cpp}
 *    ((void)((a),0))
 *
 *    ((void)((a),(b),0))
 *
 *    ((void)((a),(b),(c),0))
 * \endcode
 *
 * corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
 */
#  define AMREX_D_EXPR(a,b,c) ((void)((a),0))

/**
 * The AMREX_D_DECL(a,b,c) macro expands to a comma-separated list of
 * 1, 2, or all 3 of the arguments of the call, depending on the
 * value of AMREX_SPACEDIM.  This can be used to write portable
 * function calls that depend on 1, 2, or 3 arguments, corresponding
 * to the value of AMREX_SPACEDIM.
 *
 * For example:
 * \code{.cpp}
 *    return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
 * \endcode
 *
 * More explicitly, this macro evaluates to one of
 *
 * \code{.cpp}
 *    a
 *
 *    a,b
 *
 *    a,b,c
 * \endcode
 *
 * corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
 */
#  define AMREX_D_DECL(a,b,c) a

/**
 * \def AMREX_D_TERM
 *
 * The AMREX_D_TERM(a,b,c) macro expands to a whitespace-separated list of
 * 1, 2, or all 3 of the arguments of the call, depending on the value
 * of AMREX_SPACEDIM.  This can be used to write program logic
 * that depend on 1, 2, or 3 arguments, corresponding to
 * the value of AMREX_SPACEDIM.
 *
 * For example:
 * \code{.cpp}
 *    return AMREX_D_TERM(len[0], *len[1], *len[2]);
 * \endcode
 *
 * More explicitly, this macro evaluates to one of
 *
 * \code{.cpp}
 *    a
 *
 *    a b
 *
 *    a b c
 * \endcode
 *
 * corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
 */
#  define AMREX_D_TERM(a,b,c) a

/**
 * The AMREX_D_PICK(a,b,c) macro expands to a single result equal to the
 * 1st, 2nd, or 3rd of the arguments of the call, depending on the value
 * of AMREX_SPACEDIM.  This can be used to write program logic
 * that depend on the value of AMREX_SPACEDIM.
 *
 * For example:
 * \code{.cpp}
 *     maxsize = AMREX_D_PICK(1024, 128, 32);
 * \endcode
 *
 *  More explicitly, this macro evaluates to one of
 *
 * \code{.cpp}
 *     maxsize = 1024;
 *
 *     maxsize = 128;
 *
 *     maxsize = 32;
 * \endcode
 *
 *  corresponding to AMREX_SPACEDIM values of 1, 2, or 3.
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

/**
 * \def AMREX_1D_ONLY(a)
 * \brief Evaluate \p a only in 1-D builds; otherwise compiles to a no-op.
 *
 * \def AMREX_2D_ONLY(a)
 * \brief Evaluate \p a only in 2-D builds; otherwise compiles to a no-op.
 *
 * \def AMREX_3D_ONLY(a)
 * \brief Evaluate \p a only in 3-D builds; otherwise compiles to a no-op.
 */
#if (AMREX_SPACEDIM == 1)
#  define AMREX_1D_ONLY(a) a
#  define AMREX_2D_ONLY(a) ((void)0)
#  define AMREX_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 2)
#  define AMREX_1D_ONLY(a) ((void)0)
#  define AMREX_2D_ONLY(a) a
#  define AMREX_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 3)
#  define AMREX_1D_ONLY(a) ((void)0)
#  define AMREX_2D_ONLY(a) ((void)0)
#  define AMREX_3D_ONLY(a) a
#endif

#ifdef __cplusplus

/**
 * \def AMREX_1D_ONLY_ARGS(...)
 * \brief Expand to the variadic argument list only for 1-D builds.
 *
 * \def AMREX_2D_ONLY_ARGS(...)
 * \brief Expand to the variadic argument list only for 2-D builds.
 *
 * \def AMREX_3D_ONLY_ARGS(...)
 * \brief Expand to the variadic argument list only for 3-D builds.
 */
#if (AMREX_SPACEDIM == 1)
#define AMREX_1D_ONLY_ARGS(...) __VA_ARGS__ ,
#else
#define AMREX_1D_ONLY_ARGS(...)
#endif

#if (AMREX_SPACEDIM == 2)
#define AMREX_2D_ONLY_ARGS(...) __VA_ARGS__ ,
#else
#define AMREX_2D_ONLY_ARGS(...)
#endif

#if (AMREX_SPACEDIM == 3)
#define AMREX_3D_ONLY_ARGS(...) __VA_ARGS__ ,
#else
#define AMREX_3D_ONLY_ARGS(...)
#endif

#endif

#endif /*AMREX_SPACE_H*/
