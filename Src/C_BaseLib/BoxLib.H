#ifndef BL_BOXLIB_H
#define BL_BOXLIB_H

//
// $Id: BoxLib.H,v 1.13 2001-07-20 17:01:44 car Exp $
//

//
//@Man:
//@Memo: Error-reporting Functions and Version Number
/*@Doc:

  This class consists entirely of static error-reporting functions that
  are used throughout the BoxLib implementation together with a version
  string.  In the best of all possible worlds we'd use a namespace, but
  we aren't there, yet :-)  Note that all the output functions add a
  terminating exclamation mark, so there is no need to add any punctuation
  to sentences that get output explicitly.
*/

namespace BoxLib
{
    //
    //@ManDoc: Print out message to cerr and exit via BoxLib::Abort().
    //
    void Error (const char * msg = 0);
    //
    //@ManDoc: Print out warning message to cerr.
    //
    void Warning (const char * msg);
    //
    //@ManDoc: Print out message to cerr and exit via abort().
    //
    void Abort (const char * msg = 0);

    /*@ManDoc: Prints assertion failed messages to cerr and exits
               via abort().  Intended for use by the BL\_ASSERT() macro
               in <BLassert.H>.
    */
    void Assert (const char* EX,
                 const char* file,
                 int         line);

    /*@ManDoc: Prints out an out-of-memory message and aborts.  It is
             called by various BoxLib routines when a call to new fails.

             Called as BoxLib::OutOfMemory(\_\_FILE\_\_, \_\_LINE\_\_);
    */
    void OutOfMemory (const char* file,
                      int         line);
#ifdef NDEBUG
    const bool NDebug = true;
#else
    const bool NDebug = false;
#endif
}

#endif /*BL_BOXLIB_H*/
