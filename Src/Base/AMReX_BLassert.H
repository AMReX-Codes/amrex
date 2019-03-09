
#ifndef AMREX_BLASSERT_H_
#define AMREX_BLASSERT_H_

#include <AMReX.H>

#ifdef AMREX_ASSERT
#undef AMREX_ASSERT
#endif

#ifdef BL_ASSERT
#undef BL_ASSERT
#endif

#if !defined(AMREX_DEBUG) && !defined(AMREX_USE_ASSERTION)

/*
  The macro AMREX_ASSERT(EX) is a debugging macro, used to test the
  internal consistency of AMReX, that may also be useful in user
  code.  The argument EX is a C++ expression that you expect to
  always evaluate to true.

  When code is compiled with the DEBUG macro set to TRUE, the
  AMREX_ASSERT(EX) macro expands to a boolean test of the expression EX.
  If that expression evaluates to true, nothing happens, but if the
  expression evaluates to false, a message is output detailing the
  file and line number of the AMREX_ASSERT(EX) statement, as well as the
  literal expression EX itself, and then exits via abort() using
  amrex::Assert().  The idea is that if the assertion fails, something
  has gone terribly wrong somewhere.

  If the DEBUG macro is not set to TRUE, the AMREX_ASSERT(EX) call becomes
  a null statement in the code, which will be eliminated by the code
  optimizer.
*/
#define AMREX_ASSERT_WITH_MESSAGE(EX,MSG) ((void)0)
#define AMREX_ASSERT(EX) ((void)0)
#define    BL_ASSERT(EX) ((void)0)

#else

#define AMREX_ASSERT_WITH_MESSAGE(EX,MSG) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__ , # MSG)
#define AMREX_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)
#define    BL_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)

#endif

#define AMREX_ALWAYS_ASSERT_WITH_MESSAGE(EX,MSG) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__ , # MSG)
#define AMREX_ALWAYS_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)


#if !defined(AMREX_USE_GPU) || !defined(AMREX_DEBUG) && !defined(AMREX_USE_ASSERTION)
#define AMREX_GPU_ASSERT(EX) ((void)0)
#else
#define AMREX_GPU_ASSERT(EX) (EX)?((void)0):amrex::Assert( # EX , __FILE__, __LINE__)
#endif

#endif /*BL_BL_ASSERT_H*/
