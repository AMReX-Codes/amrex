
#ifndef BL_BOXLIB_H
#define BL_BOXLIB_H

#include <ccse-mpi.H>

#include <string>
#include <stdio.h>

//
// Error-reporting Functions and Version Number
//
// This class consists entirely of static error-reporting functions that
// are used throughout the BoxLib implementation together with a version
// string.  In the best of all possible worlds we'd use a namespace, but
// we aren't there, yet :-)  Note that all the output functions add a
// terminating exclamation mark, so there is no need to add any punctuation
// to sentences that get output explicitly.
//

namespace BoxLib
{
    void Initialize (int& argc, char**& argv, bool build_parm_parse=true, MPI_Comm mpi_comm = MPI_COMM_WORLD);

    void Finalize(bool finalize_parallel = true);
    //
    // A pointer to a function taking and returning void.
    //
    // Functions of this type, which include static class functions of
    // this type, are what we pass to ExecOnFinalize().
    //
    typedef void (*PTR_TO_VOID_FUNC)(void);
    //
    // We maintain a stack of functions that need to be called in Finalize().
    // The functions are called in LIFO order.  The idea here is to allow
    // classes to clean up any "global" state that they maintain when we're
    // exiting from BoxLib land.
    //
    void ExecOnFinalize (PTR_TO_VOID_FUNC);
    void ExecOnInitialize (PTR_TO_VOID_FUNC);
    //
    // Print out message to cerr and exit via BoxLib::Abort().
    //
    void Error (const char * msg = 0);
    //
    // Print out warning message to cerr.
    //
    void Warning (const char * msg);
    //
    // Print out message to cerr and exit via abort().
    //
    void Abort (const char * msg = 0);
    //
    // Prints assertion failed messages to cerr and exits
    // via abort().  Intended for use by the BL_ASSERT() macro
    // in <BLassert.H>.
    //
    void Assert (const char* EX,
                 const char* file,
                 int         line);
    //
    // Prints out an out-of-memory message and aborts.  It is
    //         called by various BoxLib routines when a call to new fails.
    //
    //        Called as BoxLib::OutOfMemory(__FILE__, __LINE__);
    //
    void OutOfMemory (const char* file,
                      int         line);

    //
    // This is used by BoxLib::Error(), BoxLib::Abort(), and BoxLib::Assert()
    // to ensure that when writing the message to stderr, that no additional
    // heap-based memory is allocated.
    //
    void write_to_stderr_without_buffering (const char* str);

#ifdef NDEBUG
    const bool NDebug = true;
#else
    const bool NDebug = false;
#endif

    extern std::string exename;

    extern int verbose;
}

#endif /*BL_BOXLIB_H*/
