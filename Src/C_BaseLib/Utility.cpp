//BL_COPYRIGHT_NOTICE

//
// $Id: Utility.cpp,v 1.16 1997-12-14 23:34:59 lijewski Exp $
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
get_initial_wall_clock_time ()
{
    struct timeval tp;

    if (gettimeofday(&tp, 0) != 0)
        BoxLib::Abort("get_time_of_day(): gettimeofday() failed");

    return tp.tv_sec + tp.tv_usec/1000000.0;
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
extern double BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
Utility::wsecond (double* t)
{
    struct timeval tp;

    gettimeofday(&tp,0);

    double dt = tp.tv_sec + tp.tv_usec/1000000.0 - BL_Initial_Wall_Clock_Time;

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

static
double
get_initial_wall_clock_time ()
{
    return RTC();
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
extern double BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
Utility::wsecond (double* t_)
{
    double t = RTC() - BL_Initial_Wall_Clock_Time;
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
        start = finish;

    double rr = double(finish - start)/CLOCKS_PER_SEC;

    if (r)
        *r = rr;

    return rr;
}

static
time_t
get_initial_wall_clock_time ()
{
    return ::time(0);
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
extern time_t BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
Utility::wsecond (double* r)
{
    time_t finish;

    time(&finish);

    double rr = double(finish - BL_Initial_Wall_Clock_Time);

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

void
Utility::OutOfMemory ()
{
    BoxLib::Error("Sorry, out of memory, bye ...");
}

//
// Ths following are the original random variations from
// "Random Number Generators: Good Ones Are Hard To Find",
// Stephen K. Park and Keith W. Miller, "Communications of the ACM",
// October 1988 Volume 31 Number 10, pages 1192 - 1201.
//
// This requires INT_MAX of 2**31 - 1 or larger.  If INT_MAX doesn't
// satisfy this but LONG_MAX does, replace all ints below with longs.

static int The_Random_Seed = 1;

void
Utility::InitRandom (int seed)
{
    assert(seed != 0);
    The_Random_Seed = seed;
}

double
Utility::Random ()
{
    const int A  = 16807;
    const int M  = 2147483647;  /* Mersenne prime 2^31-1 */
    const int Q  = 127773;      /* m div a               */
    const int R  = 2836;        /* m mod a               */

    int hi = The_Random_Seed / Q;

    The_Random_Seed = A * (The_Random_Seed - hi*Q) - R * hi;

    if (The_Random_Seed <  0)
	The_Random_Seed += M;

    return double(The_Random_Seed) / double(M);
}

#if 0
//
// A simple test of the random number generator.
// Report whether or not you've got a "good" random number generator.
//
int
main ()
{
    double u;
    for (int i = 1; i <= 10000; i++)
    {
        u = Utility::Random();
    }
    printf("The current value of seed is %d.", The_Random_Seed);
    printf("\nIt should be 1043618065.\n");
    return 0;
}
#endif

//
// Fortran entry point for Utility::Random().
//
#ifdef BL_FORT_USE_UPPERCASE

extern "C" void UTILRAND (Real* rn);

void
UTILRAND (Real* rn)
{
    assert(rn != 0);
    *rn = Utility::Random();
}

#else

extern "C" void utilrand_ (Real* rn);

void
utilrand_ (Real* rn)
{
    assert(rn != 0);
    *rn = Utility::Random();
}

#endif /*BL_FORT_USE_UPPERCASE*/


