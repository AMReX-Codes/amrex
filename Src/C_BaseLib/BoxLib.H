
#ifndef BL_BOXLIB_H
#define BL_BOXLIB_H

//
// $Id: BoxLib.H,v 1.9 2001-03-28 21:08:52 car Exp $
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
namespace BL_NAMESPACE
{
#endif

class BoxLib
{
protected:
    //
    // A null string used as a default argument by the following functions.
    //
    static const char* nullString;
public:
    //
    //@ManDoc: Print out message to cerr and exit via BoxLib::Abort().
    //
    static void Error (const char * msg = nullString);
    //
    //@ManDoc: Print out warning message to cerr.
    //
    static void Warning (const char * msg = nullString);
    //
    //@ManDoc: Print out message to cerr and exit via abort().
    //
    static void Abort (const char * msg = nullString);

    /*@ManDoc: Prints assertion failed messages to cerr and exits
               via abort().  Intended for use by the BL\_ASSERT() macro
               in <BLassert.H>.
    */
    static void Assert (const char* EX,
                        const char* file,
                        int         line);

    /*@ManDoc: Prints out an out-of-memory message and aborts.  It is
             called by various BoxLib routines when a call to new fails.

             Called as BoxLib::OutOfMemory(\_\_FILE\_\_, \_\_LINE\_\_);
    */
    static void OutOfMemory (const char* file,
                             int         line);

    /*@ManDoc: The BoxLib version string.
    
       Has the form:
    
         boxlib version 2.0 built Jun 25 1996 at 14:52:36
   */
    static const char * const version;

};

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_BOXLIB_H*/
