//BL_COPYRIGHT_NOTICE

//
// $Id: Utility.cpp,v 1.11 1997-11-25 19:09:42 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <cstring>
#include <cctype>
#else
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <REAL.H>
#include <Misc.H>
#include <BoxLib.H>
#include <Utility.H>

#ifdef WIN32
#include <direct.h>
#endif

#if !defined(BL_ARCH_CRAY) && !defined(WIN32)

#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/param.h>
#include <unistd.h>

//
// This doesn't seem to be defined on SunOS when using g++.
//
#if defined(__GNUG__) && defined(__sun) && defined(BL_SunOS)
extern "C" int gettimeofday (struct timeval*, struct timezone*);
#endif

double
Utility::second (double* t)
{
    struct tms buffer;

    times(&buffer);
    static long CyclesPerSecond = 0;
    if (CyclesPerSecond == 0)
    {
#if defined(_SC_CLK_TCK)
        CyclesPerSecond = sysconf(_SC_CLK_TCK);
        if (CyclesPerSecond == -1)
            BoxLib::Error("second(double*): sysconf() failed");
#elif defined(HZ)
        CyclesPerSecond = HZ;
#else
        CyclesPerSecond = 100;
        BoxLib::Warning("second(): sysconf(): default value of 100 for hz, worry about timings");
#endif
    }
    double dt = (buffer.tms_utime + buffer.tms_stime)/(1.0*CyclesPerSecond);
    if (t != 0)
        *t = dt;
    return dt;
}

static
double
get_time_of_day ()
{
    struct timeval tp;

    if (gettimeofday(&tp, 0) != 0)
        BoxLib::Abort("wsecond_int: gettimeofday() failed");

    return tp.tv_sec + tp.tv_usec/1000000.0;
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
static double Initial_Wall_Clock_Time = get_time_of_day();

double
Utility::wsecond (double* t)
{
    struct timeval tp;

    gettimeofday(&tp,0);

    double dt = tp.tv_sec + tp.tv_usec/1000000.0 - Initial_Wall_Clock_Time;

    if (t != 0)
        *t = dt;

    return dt;
}

#elif defined(BL_ARCH_CRAY)
extern "C" double SECOND();
extern "C" double RTC();

double
Utility::second (double* t_)
{
    double t = SECOND();
    if (t_)
        *t_ = t;
    return t;
}

double
Utility::wsecond (double* t_)
{
    static double epoch = -1.0;

    if (epoch < 0.0)
    {
        epoch = RTC();
    }

    double t = RTC() - epoch;

    if (t_)
        *t_ = t;

    return t;
}

#else

#include <time.h>

double
Utility::second (double* r)
{
    static clock_t start = -1;

    clock_t finish = clock();

    if (start == -1)
    {
	start = finish;
    }

    double rr = double(finish - start)/CLOCKS_PER_SEC;

    if (r)
        *r = rr;

    return rr;
}
double
Utility::wsecond (double* r)
{
    static time_t start = -1;

    time_t finish;

    time(&finish);

    if (start == -1)
    {
	start = finish;
    }

    double rr = double(finish - start);

    if (r)
        *r = rr;

    return rr;
}

#endif /*!defined(BL_ARCH_CRAY)*/

//
// Return true if argument is a non-zero length string of digits.
//

bool
Utility::is_integer (const char* str)
{
    if (str == 0)
        return false;

    int len = strlen(str);

    if (len == 0)
        return false;
    else
    {
        for (int i = 0; i < len; i++)
            if (!isdigit(str[i]))
                return false;
        return true;
    }
}

//
// Creates the specified directories.  `path' may be either a full pathname
// or a relative pathname.  It will create all the directories in the
// pathname, if they don't already exist, so that on successful return the
// pathname refers to an existing directory.  Returns true or false
// depending upon whether or not all it was successful.  Also returns
// true if `path' is NULL.  `mode' is the mode passed to mkdir() for
// any directories that must be created.
//
// For example, if it is passed the string "/a/b/c/d/e/f/g", it
// will return successfully when all the directories in the pathname
// exist; i.e. when the full pathname is a valid directory.
//

bool
#ifdef WIN32
Utility::CreateDirectory (const aString& path,
			  int)
#else
Utility::CreateDirectory (const aString& path,
                          mode_t         mode)
#endif
{
    if (path.length() == 0 || path == "/")
        return true;

    if (strchr(path.c_str(), '/') == 0)
    {
        //
        // No slashes in the path.
        //
#ifdef WIN32
	return _mkdir(path.c_str()) < 0 && errno != EACCES ? false : true;
#else
        return mkdir(path.c_str(),mode) < 0 && errno != EEXIST ? false : true;
#endif
    }
    else
    {
        //
        // Make copy of the directory pathname so we can write to it.
        //
        char* dir = new char[path.length() + 1];
        if (dir == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        (void) strcpy(dir, path.c_str());

        char* slash = strchr(dir, '/');

        if (dir[0] == '/')
        {
            //
            // Got a full pathname.
            //
            do
            {
                if (*(slash+1) == 0)
                    break;
                if ((slash = strchr(slash+1, '/')) != 0)
                    *slash = 0;
#ifdef WIN32
		if (_mkdir(dir) < 0 && errno != EACCES )
#else
                if (mkdir(dir, mode) < 0 && errno != EEXIST)
#endif
                    return false;
                if (slash)
                    *slash = '/';
            } while (slash);
        }
        else
        {
            //
            // Got a relative pathname.
            //
            do
            {
                *slash = 0;
#ifdef WIN32
		if (_mkdir(dir) < 0 && errno != EACCES )
#else
                if (mkdir(dir, mode) < 0 && errno != EEXIST)
#endif
                    return false;
                *slash = '/';
            } while ((slash = strchr(slash+1, '/')) != 0);

#ifdef WIN32
	    if (_mkdir(dir) < 0 && errno != EACCES)
#else
            if (mkdir(dir, mode) < 0 && errno != EEXIST)
#endif
                return false;
        }

        delete [] dir;

        return true;
    }
}

void
Utility::CreateDirectoryFailed (const aString& dir)
{
    aString msg("Couldn't create directory: ");
    msg += dir;
    BoxLib::Error(msg.c_str());
}

void
Utility::FileOpenFailed (const aString& file)
{
    aString msg("Couldn't open file: ");
    msg += file;
    BoxLib::Error(msg.c_str());
}

void
Utility::UnlinkFile (const aString& file)
{
    unlink(file.c_str());
}
