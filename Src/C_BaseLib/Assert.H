//BL_COPYRIGHT_NOTICE

#ifndef BL_ASSERT_H
#define BL_ASSERT_H

//
// $Id: Assert.H,v 1.3 1999-05-10 17:18:43 car Exp $
//

#include <BoxLib.H>

#ifdef assert
#undef assert
#endif

#if defined(NDEBUG)

/*@ManDoc:

  The macro assert(EX) is a debugging macro, used to test the
  internal consistency of BoxLib, that may also be useful in user
  code.  The argument EX is a C++ expression that you expect to
  always evaluate to true.

  When code is compiled with the DEBUG macro set to TRUE, the
  assert(EX) macro expands to a boolean test of the expression EX.
  If that expression evaluates to true, nothing happens, but if the
  expression evaluates to false, a message is output detailing the
  file and line number of the assert(EX) statement, as well as the
  literal expression EX itself, and then exits via abort() using
  BoxLib::Assert().  The idea is that if the assertion fails, something
  has gone terribly wrong somewhere.

  If the DEBUG macro is not set to TRUE, the assert(EX) call becomes
  a null statement in the code, which will be eliminated by the code
  optimizer.
*/
#define assert(EX) ((void)0)
#else
#if defined(__sgi) || defined(__hpux)
//
// The problem here is that some brain-dead compilers cannot deal with a
// void argument to the ternary operator.  Currently SGI's OCC and CC on
// HP-UX 9.0.1 have this problem.
//
#define assert(EX) if ((EX)) ; else BoxLib::Assert( # EX , __FILE__, __LINE__)
#else
#define assert(EX) (EX)?((void)0):BoxLib::Assert( # EX , __FILE__, __LINE__)
#endif
#endif

#endif /*BL_ASSERT_H*/
