#ifndef BL_BOXLIB_H
#define BL_BOXLIB_H

//
// $Id: BoxLib.H,v 1.11 2001-07-18 23:07:06 lijewski Exp $
//

//
//@Man:
//@Memo: Error-reporting Functions and Version Number
/*@Doc:

  This class consists entirely of static error-reporting functions that
  are used throughout the BoxLib implementation together with a version
  string.  In the best of all possible worlds we'd use a namespace, but
  we aren't there, yet :-)  Note that all the output functions add a
  terminating exclamation mark, so there is no need to add any punction
  to sentences that get output explicitly.
*/

#ifdef BL_NAMESPACE
namespace BoxLib
#define BL_STATIC
#else
#define BL_STATIC static
struct BoxLib
#endif
{
    //
    //@ManDoc: Print out message to cerr and exit via BoxLib::Abort().
    //
    BL_STATIC void Error (const char * msg = 0);
    //
    //@ManDoc: Print out warning message to cerr.
    //
    BL_STATIC void Warning (const char * msg);
    //
    //@ManDoc: Print out message to cerr and exit via abort().
    //
    BL_STATIC void Abort (const char * msg = 0);

    /*@ManDoc: Prints assertion failed messages to cerr and exits
               via abort().  Intended for use by the BL\_ASSERT() macro
               in <BLassert.H>.
    */
    BL_STATIC void Assert (const char* EX,
                           const char* file,
                           int         line);

    /*@ManDoc: Prints out an out-of-memory message and aborts.  It is
             called by various BoxLib routines when a call to new fails.

             Called as BoxLib::OutOfMemory(\_\_FILE\_\_, \_\_LINE\_\_);
    */
    BL_STATIC void OutOfMemory (const char* file,
                                int         line);


#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BoxLib;
#endif
#else
};
#endif

#undef BL_STATIC

#endif /*BL_BOXLIB_H*/
