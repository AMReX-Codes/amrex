//BL_COPYRIGHT_NOTICE

//
// $Id: BoxLib.cpp,v 1.2 1997-09-18 20:12:45 lijewski Exp $
//

#ifdef _MSC_VER
#include <strstrea.h>
#else
#include <strstream.h>
#endif

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#else
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#include <BoxLib.H>
#include <BLVERSION.H>

//
// The definition of our NULL string used as default argument.
//
const char* BoxLib::nullString = "";

#define bl_str(s)  # s
#define bl_xstr(s) bl_str(s)

//
// The definition of our version string.
//
const char * const BoxLib::version =

"boxlib version "
bl_xstr(BL_VERSION_MAJOR)
"."
bl_xstr(BL_VERSION_MINOR)
" built "
__DATE__
" at "
__TIME__;

#undef bl_str
#undef bl_xstr

//
// This is used by BoxLib::Error(), BoxLib::Abort(), and BoxLib::Assert()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//

static
void
write_to_stderr_without_buffering (const char* str)
{
    if (str)
    {
        //
        // Flush all buffers.
        //
        fflush(NULL);
        //
        // Add some `!'s and a newline to the string.
        //
        const char * const end = " !!!\n";
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

void
BoxLib::Error (const char* msg)
{
    write_to_stderr_without_buffering(msg);
    abort();
}

void
BoxLib::Abort (const char* msg)
{
    write_to_stderr_without_buffering(msg);
    abort();
}

void
BoxLib::Warning (const char* msg)
{
    if (msg)
        cerr << msg << '!' << endl;
}

//
// A pre-allocated buffer used by BoxLib::Assert().
//
const int DIMENSION = 512;
static char buf[DIMENSION];

void
BoxLib::Assert (const char* EX,
                const char* file,
                int         line)
{
    //
    // Why declare this static?  Well, the expectation is that by declaring
    // it static, the space for it'll be pre-allocated, so that when this
    // function is called, only the appropriate constructor will need to
    // be called.  This way BoxLib::Assert() should work fine, even when
    // complaining about running out of memory, or some such nasty situation.
    //
    static ostrstream os(buf, DIMENSION);

    os << "Assertion `"
       << EX
       << "' failed, "
       << "file \""
       << file
       << "\", "
       << "line "
       << line
       << ends;

    write_to_stderr_without_buffering(buf);

    abort();
}

void
BoxLib::OutOfMemory (const char* file,
                     int         line)
{
    BoxLib::Assert("operator new", file, line);
}
