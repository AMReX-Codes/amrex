#ifndef BL_UTILITY_H
#define BL_UTILITY_H

//
// $Id: Utility.H,v 1.34 2001-07-17 23:02:28 lijewski Exp $
//

#include <iostream>

#ifndef WIN32
#include <sys/types.h>
#include <sys/wait.h>
#endif

#include <aString.H>
#include <REAL.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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
               relative to program startup.  If t != 0, it also
               stores the result in t.
    */
    static double wsecond (double* t = 0);
    //
    //@ManDoc: Reset start of Wall Clock Time for wsecond() to NOW.
    //
    static void ResetWallClockTime ();
    //
    //@ManDoc: Return true if argument is a non-zero length string of digits.
    //
    static bool is_integer (const char* str);

#ifdef BL_DOES_USE_BSP_NOW
    //
    // The `Op' functions aren't documented as they're used only in BSP-specific code.
    //
    static void OpBoolAnd (bool& result,
                           bool& arg1,
                           bool& arg2,
                           int*  nbytes)
    {
        result = arg1 && arg2;
    }
    static void OpBoolOr  (bool& result,
                           bool& arg1,
                           bool& arg2,
                           int*  nbytes)
    {
        result = arg1 || arg2;
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
                           int*  nbytes)
    {
        result = arg1 && arg2;
    }
#endif
    //
    //@ManDoc: Returns rootNNNN where NNNN == num.
    //
    static aString Concatenate (const aString& root,
                                int            num);

#ifdef WIN32
    static bool UtilCreateDirectory (const aString& path,
                                     int            mode);
#else
    /*@ManDoc:
      Creates the specified directories.  `path' may be either a full pathname
      or a relative pathname.  It will create all the directories in the
      pathname, if they don't already exist, so that on successful return the
      pathname refers to an existing directory.  Returns true or false
      depending upon whether or not all it was successful.  Also returns
      if `path' is NULL or "/".  `mode' is the mode passed to mkdir() for
      any directories that must be created.
    
      For example, if it is passed the string "/a/b/c/d/e/f/g", it
      will return successfully when all the directories in the pathname
      exist; i.e. when the full pathname is a valid directory.

      In a Windows environment, the path separator is a '\', so that if using
      the example given above you must pass the string  
      "\\a\\b\\c\\d\\e\\f\\g"  (Note that you must escape the backslash in a
      character string), 
    */
    static bool UtilCreateDirectory (const aString& path,
                                     mode_t         mode);
#endif
    //
    //@ManDoc: Output a message and abort when couldn't create the directory.
    //
    static void CreateDirectoryFailed (const aString& dir);
    //
    //@ManDoc: Output a message and abort when couldn't open the file.
    //
    static void FileOpenFailed (const aString& dir);
    //
    //@ManDoc: Attempt to unlink the file.  Ignore any errors.
    //
    static void UnlinkFile (const aString& file);

    /*@ManDoc: Aborts after printing message indicating out-of-memory;
               i.e. operator new has failed. This is the "supported"
               set\_new\_handler() function for BoxLib applications.
     */
    static void OutOfMemory ();

    /*@ManDoc:
      Mersenne Twister pseudo-random number generator.

      Generates one pseudorandom real number (double) which is
      uniformly distributed on [0,1]-interval for each call.

      Accepts any 32-bit integer as a seed -- uses 4357 as the default.

      Has a period of 2**19937.

      Mersenne Twister Home Page: http://www.math.keio.ac.jp/matumoto/emt.html

      There is also an entry point for Fortran callable as:
    
      REAL\_T rn
      call blutilrand(rn)
    */
    static double Random ();

    /*@ManDoc: Set the seed of the random number generator.
    */
    static void InitRandom (unsigned long seed);
#ifndef WIN32
    /*@ManDoc: Try to execute a shell command in a subshell.
               User is responsible for calling waitpid() on the pid
	       returned from this function.
    */
    static pid_t Execute (const char* cmd);
#endif
};

//
// I'm going to document right here all the BL macros that aren't documented
// anywhere else.  Note that all these #ifdef ... #endif blocks are necessary
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
#endif /*BL_FORT_USE_UNDERSCORE*/

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

#ifdef BL_FORT_USE_LOWERCASE
#undef BL_FORT_USE_LOWERCASE
/*@ManDoc:
  The macro BL\_FORT\_USE\_LOWERCASE indicates that C++ code should call
  Fortran routines using lowercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.  

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyx.
*/
#define BL_FORT_USE_LOWERCASE 1
#endif /*BL_FORT_USE_LOWERCASE*/

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

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif


#endif /*BL_UTILITY_H*/
