
//
// $Id: Utility.cpp,v 1.42 2000-10-03 20:26:49 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdio>
#else
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/wait.h>
#endif
#include <errno.h>

#include <REAL.H>
#include <Misc.H>
#include <BoxLib.H>
#include <Utility.H>

#ifdef WIN32
#include <direct.h>
#define mkdir(a,b) _mkdir((a))
const char* path_sep_str = "\\";
#else
const char* path_sep_str = "/";
#endif

#ifdef BL_T3E
#include <malloc.h>
#endif

#if !defined(BL_ARCH_CRAY) && !defined(WIN32) && !defined(BL_T3E)

#include <sys/types.h>
#include <sys/times.h>
#ifdef BL_AIX
#undef _XOPEN_SOURCE_EXTENDED
#define _XOPEN_SOURCE_EXTENDED
#endif
#include <sys/time.h>
#ifdef BL_AIX
#undef _XOPEN_SOURCE_EXTENDED
#endif
#include <sys/param.h>
#include <unistd.h>

//
// This doesn't seem to be defined on SunOS when using g++.
//
#if defined(__GNUG__) && defined(__sun) && defined(BL_SunOS)
extern "C" int gettimeofday (struct timeval*, struct timezone*);
#endif

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
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
double BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

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

#include <unistd.h>

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
double BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
Utility::wsecond (double* t_)
{
    double t = RTC() - BL_Initial_Wall_Clock_Time;
    if (t_)
        *t_ = t;
    return t;
}

#elif defined(BL_T3E)

//#include <intrinsics.h>
#include <unistd.h>
extern "C" long _rtc();

static double BL_Clock_Rate;
extern "C"
{
long IRTC_RATE();
long _irt();
}

static
long
get_initial_wall_clock_time ()
{
    BL_Clock_Rate = IRTC_RATE();
    return _rtc();
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
long BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

//
// NOTE: this is returning wall clock time, instead of cpu time.  But on
// the T3E, there is no difference (currently).  If we call second() instead,
// we may be higher overhead.  Think about this one.
//
double
Utility::second (double* t_)
{
    double t = (_rtc() - BL_Initial_Wall_Clock_Time)/BL_Clock_Rate;
    if (t_)
        *t_ = t;
    return t;
}

double
Utility::wsecond (double* t_)
{
    double t = (_rtc() - BL_Initial_Wall_Clock_Time)/BL_Clock_Rate;
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
time_t BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

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

#endif /*!defined(BL_ARCH_CRAY) && !defined(WIN32) && !defined(BL_T3E)*/

void
Utility::ResetWallClockTime ()
{
    BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();
}

//
// Return true if argument is a non-zero length string of digits.
//

bool
Utility::is_integer (const char* str)
{
    int len = 0;

    if (str == 0 || (len = strlen(str)) == 0)
        return false;

    for (int i = 0; i < len; i++)
        if (!isdigit(str[i]))
            return false;

    return true;
}

aString
Utility::Concatenate (const aString& root,
                      int            num)
{
    aString result = root;
    char buf[sizeof(int) + 1];
    sprintf(buf, "%04d", num);
    result += buf;
    return result;
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
Utility::UtilCreateDirectory (const aString& path,
                              int)
#else
Utility::UtilCreateDirectory (const aString& path,
                              mode_t         mode)
#endif
{
    if (path.length() == 0 || path == path_sep_str)
        return true;

    if (strchr(path.c_str(), *path_sep_str) == 0)
    {
        //
        // No slashes in the path.
        //
        return mkdir(path.c_str(),mode) < 0 && errno != EEXIST ? false : true;
    }
    else
    {
        //
        // Make copy of the directory pathname so we can write to it.
        //
        char* dir = new char[path.length() + 1];
        (void) strcpy(dir, path.c_str());

        char* slash = strchr(dir, *path_sep_str);

        if (dir[0] == *path_sep_str)
        {
            //
            // Got a full pathname.
            //
            do
            {
                if (*(slash+1) == 0)
                    break;
                if ((slash = strchr(slash+1, *path_sep_str)) != 0)
                    *slash = 0;
                if (mkdir(dir, mode) < 0 && errno != EEXIST)
                    return false;
                if (slash)
                    *slash = *path_sep_str;
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
                if (mkdir(dir, mode) < 0 && errno != EEXIST)
                    return false;
                *slash = *path_sep_str;
            } while ((slash = strchr(slash+1, *path_sep_str)) != 0);

            if (mkdir(dir, mode) < 0 && errno != EEXIST)
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
#ifdef BL_T3E
    malloc_stats(0);
#endif
    BoxLib::Error("Sorry, out of memory, bye ...");
}

/*
** Mersenne Twister pseudo-random number generator.
**
** Generates one pseudorandom real number (double) which is
** uniformly distributed on [0,1]-interval for each call.
**
** Uses any 32-bit integer except for 0 for a seed (with 4357 as the default).
**
** Has a period of 2**19937.
**
** Coded by Takuji Nishimura, considering the suggestions by
** Topher Cooper and Marc Rieffel in July-Aug. 1997.
**
** Mersenne Twister Home Page: http://www.math.keio.ac.jp/matumoto/emt.html
**
** Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
*/

#define RAND_N 624
#define RAND_M 397

static unsigned long mt[RAND_N];
static int mti = RAND_N+1;

void
Utility::InitRandom (unsigned long seed)
{
    /*
    ** Setting initial seeds using the generator Line 25 of Table 1 in
    ** KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102.
    */
    int i;

    for (i = 0; i < RAND_N; i++)
    {
         mt[i] = seed & 0xffff0000;
         seed = 69069 * seed + 1;
         mt[i] |= (seed & 0xffff0000) >> 16;
         seed = 69069 * seed + 1;
    }
    mti = RAND_N;
}

double
Utility::Random ()
{
    static unsigned long mag01[2] = { 0x0, 0x9908b0df };

    unsigned long y;

    if (mti >= RAND_N)
    {
        /*
        ** Generate RAND_N words at one time.
        */
        int kk;

        if (mti == RAND_N+1)
            /*
            ** Use the default initial seed.
            */
            Utility::InitRandom(4357);

        for (kk = 0; kk < RAND_N-RAND_M; kk++)
        {
            y = (mt[kk]&0x80000000)|(mt[kk+1]&0x7fffffff);
            mt[kk] = mt[kk+RAND_M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for ( ; kk < RAND_N-1; kk++)
        {
            y = (mt[kk]&0x80000000)|(mt[kk+1]&0x7fffffff);
            mt[kk] = mt[kk+(RAND_M-RAND_N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[RAND_N-1]&0x80000000)|(mt[0]&0x7fffffff);
        mt[RAND_N-1] = mt[RAND_M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return ( (double)y / (unsigned long)0xffffffff );
}

#undef RAND_N
#undef RAND_M

//
// Fortran entry point for Utility::Random().
//
#if defined(BL_FORT_USE_UPPERCASE)

extern "C" void BLUTILRAND (Real* rn);

void
BLUTILRAND (Real* rn)
{
    BL_ASSERT(rn != 0);
    *rn = Utility::Random();
}

#elif defined(BL_FORT_USE_LOWERCASE)
extern "C" void blutilrand (Real* rn);

void
blutilrand (Real* rn)
{
    BL_ASSERT(rn != 0);
    *rn = Utility::Random();
}

#elif defined(BL_FORT_USE_UNDERSCORE)
extern "C" void blutilrand_ (Real* rn);

void
blutilrand_ (Real* rn)
{
    BL_ASSERT(rn != 0);
    *rn = Utility::Random();
}
#endif 

#ifndef WIN32
extern "C" pid_t fork();

pid_t
Utility::Execute (const char* cmd)
{

    pid_t pid = fork();

    if (pid == 0)
    {
        system(cmd);

        exit(0);
    }

    return pid;
}
#endif

#ifdef BL_NAMESPACE
}
#endif
