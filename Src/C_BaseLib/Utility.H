/* BL_COPYRIGHT_NOTICE */

#ifndef BL_UTILITY_H
#define BL_UTILITY_H

//
// $Id: Utility.H,v 1.11 1997-11-18 18:22:24 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#ifndef WIN32
#include <sys/types.h>
#endif

#include <aString.H>
#include <Boolean.H>
#include <REAL.H>

//
//@Man:
//@Memo: Useful C++ Utility Functions
/*@Doc:

  This data-less class is a poor-person's namespace of utility functions.
  Since we can't assume the existence of namespaces, and we don't
  like global functions, we put them into a data-less class as
  static member functions.
*/

class Utility
{
public:
    /*@ManDoc: Returns the total user and system time used by
               the calling process up to the point of the call.
               If t != 0, it also stores the result in t.
    */
    static double second (double* t = 0);

    /*@ManDoc: Returns the current time as the number of microseconds
               relative to some epoch in the past.  If t != 0, it also
               stores the result in t.
    */
    static double wsecond (double* t = 0);
    //
    //@ManDoc: Return true if argument is a non-zero length string of digits.
    //
    static bool is_integer (const char* str);

    //
    // These should be templates like the following.
    // problems with void * cast to bsp.
    //
    template<class T> static void OpMin (T&   result,
                                         T&   arg1,
                                         T&   arg2,
                                         int* nbytes)
    {
        result = (arg1 < arg2 ? arg1 : arg2);
    }
    static void OpBoolAnd (bool& result,
                           bool& arg1,
                           bool& arg2,
                           int*  nbytes)
    {
        result = arg1 && arg2;
    }
    static void OpRealSum (Real& result,
                           Real& arg1,
                           Real& arg2,
                           int*  nbytes)
    {
        result = arg1 + arg2;
    }
    static void OpRealMax (Real& result,
                           Real& arg1,
                           Real& arg2,
                           int*  nbytes)
    {
        result = (arg1 > arg2 ? arg1 : arg2);
    }
    static void OpRealMin (Real& result,
                           Real& arg1,
                           Real& arg2,
                           int*  nbytes)
    {
        result = (arg1 < arg2 ? arg1 : arg2);
    }
    static void OpIntMin (int& result,
                          int& arg1,
                          int& arg2,
                          int* nbytes)
    {
        result = (arg1 < arg2 ? arg1 : arg2);
    }
    static void OpIntMax (int& result,
                          int& arg1,
                          int& arg2,
                          int* nbytes)
    {
        result = (arg1 > arg2 ? arg1 : arg2);
    }
    static void OpIntSum (int& result,
                          int& arg1,
                          int& arg2,
                          int* nbytes)
    {
        result = arg1 + arg2;
    }
    static void OpLongMin (long& result,
                           long& arg1,
                           long& arg2,
                           int*  nbytes)
    {
        result = (arg1 < arg2 ? arg1 : arg2);
    }
    static void OpLongMax (long& result,
                           long& arg1,
                           long& arg2,
                           int*  nbytes)
    {
        result = (arg1 > arg2 ? arg1 : arg2);
    }
    static void OpLongSum (long& result,
                           long& arg1,
                           long& arg2,
                           int*  nbytes)
    {
        result = arg1 + arg2;
    }
    static void OpLongAnd (long& result,
			   long& arg1,
			   long& arg2,
			   int* nbytes)
    {
        result = arg1 && arg2;
    }
    //
    // Creates the specified directories.  `path' may be either a full pathname
    // or a relative pathname.  It will create all the directories in the
    // pathname, if they don't already exist, so that on successful return the
    // pathname refers to an existing directory.  Returns true or false
    // depending upon whether or not all it was successful.  Also returns
    // if `path' is NULL or "/".  `mode' is the mode passed to mkdir() for
    // any directories that must be created.
    //
    // For example, if it is passed the string "/a/b/c/d/e/f/g", it
    // will return successfully when all the directories in the pathname
    // exist; i.e. when the full pathname is a valid directory.
    //
#ifdef WIN32
    static bool CreateDirectory (const aString& path,
				 int		mode);
#else
    static bool CreateDirectory (const aString& path,
                                 mode_t         mode);
#endif
};

//
// I'm going to document right here all the BL macros that aren't documented
// anywhere else.  Note that all there #ifdef ... #endif blocks are necessary
// to get doc++ to properly document the macros.
//

#ifdef BL_ARCH_CRAY
#undef BL_ARCH_CRAY
/*@ManDoc:
  The macro BL\_ARCH\_CRAY is defined only when compiling on a Cray
  architecture.  The make subsystem automatically sets this if it
  detects that a compilation is taking place on a Cray system.
*/
#define BL_ARCH_CRAY 1
#endif /*BL_ARCH_CRAY*/

#ifdef BL_ARCH_IEEE
#undef BL_ARCH_IEEE
/*@ManDoc:
  The macro BL\_ARCH\_IEEE is defined only when compiling on a Cray architecture.
  The make subsystem automatically sets this if it detects that a compilation
  is taking place on a non-Cray system.  That is to say, we assume all
  architectures, with the exception of Cray, are IEEE.
*/
#define BL_ARCH_IEEE 1
#endif /*BL_ARCH_IEEE*/

#ifdef BL_USE_FLOAT
#undef BL_USE_FLOAT
/*@ManDoc:
  The macro BL\_USE\_FLOAT indicates that C++ floating-point calculations should
  use "float" variables and Fortran floating-point calculations should use
  "real*4" variables.  One of BL\_USE\_FLOAT or BL\_USE\_DOUBLE must always be
  defined when compiling and using BoxLib.
*/
#define BL_USE_FLOAT 1
#endif /*BL_USE_FLOAT*/

#ifdef BL_USE_DOUBLE
#undef BL_USE_DOUBLE
/*@ManDoc:
  The macro BL\_USE\_DOUBLE indicates that C++ floating-point calculations
  should use "double" variables and Fortran floating-point calculations
  should use "real*8" variables.  One of BL\_USE\_FLOAT or BL\_USE\_DOUBLE must
  always be defined when compiling and using BoxLib.  This macro is not
  allowed on Cray architectures; i.e. only BL\_USE\_FLOAT is allowed on Cray
  architectures.
*/
#define BL_USE_DOUBLE 1
#endif /*BL_USE_DOUBLE*/

#ifdef BL_LANG_FORT
#undef BL_LANG_FORT
/*@ManDoc:
  The macro BL\_LANG\_FORT indicates that Fortran code is being compiled.
*/
#define BL_LANG_FORT 1
#endif /*BL_LANG_FORT*/

#ifdef BL_LANG_CC
#undef BL_LANG_CC
/*@ManDoc:
  The macro BL\_LANG\_CC indicates that C++ code is being compiled.
*/
#define BL_LANG_CC 1
#endif /*BL_LANG_CC*/

#ifdef BL_AUTO_INSTANTIATE
#undef BL_AUTO_INSTANTIATE
/*@ManDoc:
  The macro BL\_AUTO\_INSTANTIATE indicates that the C++ compiler is not
  able to find a template source file (.C file) given the template declaration
  file (.H file).  This macro instructs the compiler to automatically include
  the .C file, whenever it includes the .H file.
*/
#define BL_AUTO_INSTANTIATE 1
#endif /*BL_AUTO_INSTANTIATE*/


#ifdef BL_FORT_USE_UNDERSCORE
#undef BL_FORT_USE_UNDERSCORE
/*@ManDoc:
  The macro BL\_FORT\_USE\_UNDERSCORE indicates that C++ code should call
  Fortran routines by appending an underscore to the name of the Fortran
  routine.  This is set automatically by the make subsystem.  It is the
  default for all BoxLib-style programs except on Cray architectures.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyz\_.
*/
#define BL_FORT_USE_UNDERSCORE 1
#endif /*BL_LANG_CC*/

#ifdef BL_FORT_USE_UPPERCASE
#undef BL_FORT_USE_UPPERCASE
/*@ManDoc:
  The macro BL\_FORT\_USE\_UPPERCASE indicates that C++ code should call
  Fortran routines using uppercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.  Currently it
  is only set when compiling on Cray architectures.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as ABCXYZ.
*/
#define BL_FORT_USE_UPPERCASE 1
#endif /*BL_FORT_USE_UPPERCASE*/

/*@ManDoc:
  BL\_IGNORE\_MAX is a macro that expands to the literal value 100000.  It is
  defined when compiling either Fortran or C++ code; i.e. when either
  BL\_LANG\_CC or BL\_LANG\_FORT is defined.  It is used in calls to
  istream::ignore() in the BoxLib code when reading in characters from an
  istream.  We use this macro instead of the more proper INT\_MAX from
  <limits.h> since at least one compiler didn't work properly when
  istream::ignore() was passed INT\_MAX.
*/
#define BL_IGNORE_MAX 100000


#endif /*BL_UTILITY_H*/
