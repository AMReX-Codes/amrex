//BL_COPYRIGHT_NOTICE

//
// $Id: Tracer.cpp,v 1.6 2000-04-24 17:52:38 car Exp $
//
// Definition of Tracer member functions.
//

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
#include <Tracer.H>
#include <ParallelDescriptor.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

const char* Tracer::m_program = 0;

int Tracer::m_count = 0;

bool Tracer::m_trace = false;

char* Tracer::m_substr = 0;

const int IDENT = 3;  // Indent by three spaces.

Tracer::Tracer (const char* function)
{
    m_func = function;

    if (m_count++ == 0)
    {
        //
        // Only call init() the first time.
        // Also set `program'.
        // We strip off any leading '/'s from the name.
        //
        char* slash = strrchr(m_func, '/');

        m_program = slash == 0 ? m_func : ++slash;

        init();
    }

    if (m_trace && (m_substr == 0 || strstr(m_func, m_substr)))
    {
        const int DIM = 1024;

        char buf[DIM+1], *bp = buf;

        for (int i = 0; i < IDENT*m_count-1; i++)
            *bp++ = ' ';

        sprintf(bp,
                "%s (%d): entered %s\n",
                m_program,
                ParallelDescriptor::MyProc(),
                m_func);

        buf[DIM] = 0;  // Just to be really careful.

        printf(buf);
    }
}

Tracer::~Tracer ()
{
    if (m_trace && (m_substr == 0 || strstr(m_func, m_substr)))
    {
        const int DIM = 1024;

        char buf[DIM+1], *bp = buf;

        for (int i = 0; i < IDENT*m_count-1; i++)
            *bp++ = ' ';

        sprintf(bp,
                "%s (%d): exited %s\n",
                m_program,
                ParallelDescriptor::MyProc(),
                m_func);

        buf[DIM] = 0;  // Just to be really careful.

        printf(buf);
    }

    if (--m_count == 0)
    {
        delete [] m_substr;
    }
}

void Tracer::init ()
{
    char* symbol = new char[strlen(m_program) + strlen("TRACER") + 1];

    strcpy(symbol, m_program);
    strcat(symbol, "TRACER");

    char* definition = getenv(symbol);

    delete [] symbol;

    if (definition)
    {
        m_trace = true;

        if (*definition)
        {
            //
            // If it isn't the NULL string, we do function tracing.
            //
            m_substr = new char[strlen(definition)+1];

            strcpy(m_substr, definition);
        }
    }
}

#ifdef TEST_TRACER

//
// Simple test program for `Tracer'.
//
// And then run as the following to get a feel for how to use it
//
//   setenv   testTRACER;               ./tTracer
//   setenv   testTRACER dog;           ./tTracer
//   setenv   testTRACER cat;           ./tTracer
//   setenv   testTRACER worm;          ./tTracer
//
//   unsetenv testTRACER;               ./tTracer
//

static
void
dog1 ()
{
    TRACER("dog1");
}

static
void
dog2 ()
{
    TRACER("dog2");
}

static
void
cat1 ()
{
    TRACER("cat1");
    dog1();
}

static
void
cat2 ()
{
    TRACER("cat2");
    dog2();
}

int
main ()
{
    //
    // The argument to this Tracer object is taken to be the name of this
    // program.  If the environment variable "testTRACER" exists, we turn
    // tracing on.  If that environment variable has a value, we use that
    // value as an extended regular expression, and trace functions whose
    // name contains a match for the regular expression.
    //
    TRACER("test");
    cat1();
    cat2();

    return 0;
}

#endif /*TEST_TRACER*/

#ifdef BL_NAMESPACE
}
#endif
