//BL_COPYRIGHT_NOTICE

#ifndef BL_BOOLEAN_H
#define BL_BOOLEAN_H

//
// $Id: Boolean.H,v 1.3 1999-02-09 21:21:45 lijewski Exp $
//

#if !defined(BL_LANG_FORT)

#if defined(BL_HAS_NO_BOOL)

#ifdef true
#undef true
#endif

#ifdef false
#undef false
#endif

#ifdef bool
#undef bool
#endif

/*@ManDoc:
  bool is typedef'd to an int.  It is defined only for those compilers
  that do not yet support the bool type as specified in the ANSI C++ Draft
  Standard.  It is provided here as a convenience so that all C++ code can
  be written using bool as the bool type, regardless of whether or not the
  compiler truly supports the bool type.
*/
    typedef int bool;

/*@ManDoc:
  true is a constant of type bool.  It is provided here for those compilers
  that do not yet support the bool type as specified in the ANSI C++ Draft
  Standard.
*/
    const bool true = 1;

/*@ManDoc:
  false is a constant of type bool.  It is provided here for those compilers
  that do not yet support the bool type as specified in the ANSI C++ Draft
  Standard.
*/
    const bool false = 0;

#endif /*defined(BL_HAS_NO_BOOL)*/

#else

#ifdef TRUE
#undef TRUE
#endif

#ifdef FALSE
#undef FALSE
#endif

#define TRUE   1
#define FALSE  0

#endif /*!defined(BL_LANG_FORT)*/

#endif /*BL_BOOLEAN_H*/

