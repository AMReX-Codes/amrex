//BL_COPYRIGHT_NOTICE

#ifndef _BOOLEAN_H_
#define _BOOLEAN_H_

//
// $Id: Boolean.H,v 1.1 1997-09-12 18:00:02 lijewski Exp $
//

#undef TRUE
#undef FALSE
#undef True
#undef False
#undef true
#undef false
#undef bool

#if !defined(BL_LANG_FORT)
#if defined( BL_HAS_NO_BOOL )
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
#endif /* defined( BL_HAS_NO_BOOL) */

#else

#define TRUE   1
#define FALSE  0
#endif

#endif /*_BOOLEAN_H_*/

