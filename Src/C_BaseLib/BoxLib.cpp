//BL_COPYRIGHT_NOTICE

//
// $Id: BoxLib.cpp,v 1.5 1997-11-13 18:45:14 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
using std::cerr;
#else
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#endif

#include <BoxLib.H>
#include <BLVERSION.H>
#include <ParallelDescriptor.H>

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
    //
    // Flush all buffers.
    //
    fflush(NULL);

    if (str)
    {
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
#ifdef BL_USE_BSP
    ParallelDescriptor::Abort(BoxLib::nullString);
#else
    abort();
#endif
}

void
BoxLib::Abort (const char* msg)
{
    write_to_stderr_without_buffering(msg);
#ifdef BL_USE_BSP
    ParallelDescriptor::Abort(BoxLib::nullString);
#else
    abort();
#endif
}

void
BoxLib::Warning (const char* msg)
{
    if (msg)
    {
        cerr << msg << '!' << '\n';
    }
}

void
BoxLib::Assert (const char* EX,
                const char* file,
                int         line)
{
    const int DIMENSION = 1024;

    char buf[DIMENSION+1];

    sprintf(buf,
            "Assertion `%s' failed, file \"%s\", line %d",
            EX,
            file,
            line);
    //
    // Just to be a little safer :-)
    //
    buf[DIMENSION] = 0;

    write_to_stderr_without_buffering(buf);

    abort();
}

void
BoxLib::OutOfMemory (const char* file,
                     int         line)
{
    BoxLib::Assert("operator new", file, line);
}
