//
// $Id: Utility.cpp,v 1.49 2001-07-22 18:11:03 car Exp $
//

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdio>

#include <sys/stat.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/wait.h>
#endif
#include <errno.h>

#include <REAL.H>
#include <BoxLib.H>
#include <Utility.H>
#include <BLassert.H>

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

double
BoxLib::second (double* t)
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
BoxLib::wsecond (double* t)
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
BoxLib::second (double* t_)
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
BoxLib::wsecond (double* t_)
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
BoxLib::second (double* t_)
{
    double t = (_rtc() - BL_Initial_Wall_Clock_Time)/BL_Clock_Rate;
    if (t_)
        *t_ = t;
    return t;
}

double
BoxLib::wsecond (double* t_)
{
    double t = (_rtc() - BL_Initial_Wall_Clock_Time)/BL_Clock_Rate;
    if (t_)
        *t_ = t;
    return t;
}

#else

#include <time.h>

double
BoxLib::second (double* r)
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
BoxLib::wsecond (double* r)
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
BoxLib::ResetWallClockTime ()
{
    BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();
}

//
// Return true if argument is a non-zero length string of digits.
//

bool
BoxLib::is_integer (const char* str)
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
BoxLib::Concatenate (const aString& root,
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
BoxLib::UtilCreateDirectory (const aString& path,
                             int)
#else
BoxLib::UtilCreateDirectory (const aString& path,
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
BoxLib::CreateDirectoryFailed (const aString& dir)
{
    aString msg("Couldn't create directory: ");
    msg += dir;
    BoxLib::Error(msg.c_str());
}

void
BoxLib::FileOpenFailed (const aString& file)
{
    aString msg("Couldn't open file: ");
    msg += file;
    BoxLib::Error(msg.c_str());
}

void
BoxLib::UnlinkFile (const aString& file)
{
    unlink(file.c_str());
}

void
BoxLib::OutOfMemory ()
{
#ifdef BL_T3E
    malloc_stats(0);
#endif
    BoxLib::Error("Sorry, out of memory, bye ...");
}

#ifndef WIN32
extern "C" pid_t fork();

pid_t
BoxLib::Execute (const char* cmd)
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


//
// Encapsulates Time
//

namespace
{
const long billion = 1000000000L;
}

BoxLib::Time::Time()
{
    tv_sec = 0;
    tv_nsec = 0;
}

BoxLib::Time::Time(long s, long n)
{
    BL_ASSERT(s >= 0);
    BL_ASSERT(n >= 0);
    BL_ASSERT(n < billion);
    tv_sec = s;
    tv_nsec = n;
    normalize();
}

BoxLib::Time::Time(double d)
{
    tv_sec = long(d);
    tv_nsec = long((d-tv_sec)*billion);
    normalize();
}

double
BoxLib::Time::as_double() const
{
    return tv_sec + tv_nsec/double(billion);
}

long
BoxLib::Time::as_long() const
{
    return tv_sec + tv_nsec/billion;
}

BoxLib::Time&
BoxLib::Time::operator+=(const Time& r)
{
    tv_sec += r.tv_sec;
    tv_nsec += r.tv_nsec;
    normalize();
    return *this;
}

BoxLib::Time
BoxLib::Time::operator+(const Time& r) const
{
    Time result(*this);
    return result+=r;
}

void
BoxLib::Time::normalize()
{
    if ( tv_nsec > billion )
    {
	tv_nsec -= billion;
	tv_sec += 1;
    }
}

BoxLib::Time
BoxLib::Time::get_time()
{
    return Time(BoxLib::wsecond());
}


//
// BoxLib Interface to Mersenne Twistor
//

/* A C-program for MT19937: Real number version (1999/10/28)    */
/*   genrand() generates one pseudorandom real number (double)  */
/* which is uniformly distributed on [0,1]-interval, for each   */
/* call. sgenrand(seed) sets initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be       */
/* called once. (seed is any 32-bit integer.)                   */
/* Integer generator is obtained by modifying two lines.        */
/*   Coded by Takuji Nishimura, considering the suggestions by  */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.            */

/* This library is free software under the Artistic license:       */
/* see the file COPYING distributed together with this code.       */
/* For the verification of the code, its output sequence file      */
/* mt19937-1.out is attached (2001/4/2)                           */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp                                        */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

namespace
{
// Period parameters
// const int N = 624;
const int M = 397;
const unsigned long MATRIX_A   = 0x9908B0DFUL; // constant vector a
const unsigned long UPPER_MASK = 0x80000000UL; // most significant w-r bits
const unsigned long LOWER_MASK = 0x7FFFFFFFUL; // least significant r bits

// Tempering parameters
const unsigned long TEMPERING_MASK_B = 0x9D2C5680UL;
const unsigned long TEMPERING_MASK_C = 0xEFC60000UL;

inline unsigned long TEMPERING_SHIFT_U(unsigned long y) { return y >> 11L; }
inline unsigned long TEMPERING_SHIFT_S(unsigned long y) { return y << 7L ; }
inline unsigned long TEMPERING_SHIFT_T(unsigned long y) { return y << 15L; }
inline unsigned long TEMPERING_SHIFT_L(unsigned long y) { return y >> 18L; }
}

// initializing the array with a NONZERO seed
void
BoxLib::mt19937::sgenrand(unsigned long seed)
{
    for ( int i = 0; i < N; ++i )
    {
	mt[i] = seed & 0xFFFF0000UL;
	seed  = 69069U * seed + 1;
	mt[i] |= (seed& 0xFFFF0000UL) >> 16;
	seed = 69069U*seed + 1;
    }
    mti = N;
}

void
BoxLib::mt19937::reload()
{
    unsigned long y;
    int kk;
    // mag01[x] = x * MATRIX_A  for x=0,1
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    for ( kk=0; kk<N-M; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+M] ^ (y >> 1L) ^ mag01[y & 0x1];
    }
    for ( ; kk<N-1; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+(M-N)] ^ (y >> 1L) ^ mag01[y & 0x1];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1L) ^ mag01[y & 0x1];

    mti = 0;
}

unsigned long
BoxLib::mt19937::igenrand()
{
    // generate N words at one time
    if ( mti >= N ) reload();

    unsigned long y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y;
}

BoxLib::mt19937::mt19937(unsigned long seed)
    : init_seed(seed), mti(N+1)
{
    sgenrand(seed);
}

void
BoxLib::mt19937::rewind()
{
    sgenrand(init_seed);
}

double
BoxLib::mt19937::d1_value()
{
    return double(igenrand())/0xFFFFFFFFUL;
}

double
BoxLib::mt19937::d_value()
{
    const double zzz = double(0x80000000UL)*2;
    return double(igenrand())/zzz;
}

long
BoxLib::mt19937::l_value()
{
    return igenrand()&0x7FFFFFFFUL;
}

unsigned long
BoxLib::mt19937::u_value()
{
    return igenrand();
}

namespace
{
    BoxLib::mt19937 the_generator;
}

void
BoxLib::InitRandom (unsigned long seed)
{
    the_generator = mt19937(seed);
}

double
BoxLib::Random ()
{
    return the_generator.d_value();
}

//
// Fortran entry point for BoxLib::Random().
//

#if defined(BL_FORT_USE_UPPERCASE)
#define FORT_BLUTILRAND BLUTILRAND
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_BLUTILRAND blutilrand
#else defined(BL_FORT_USE_UNDERSCORE)
#define FORT_BLUTILRAND blutilrand_
#endif

extern "C" void FORT_BLUTILRAND (Real* rn);

void
FORT_BLUTILRAND (Real* rn)
{
    BL_ASSERT(rn != 0);
    *rn = BoxLib::Random();
}


//
// Sugar for parsing IO
//

std::istream&
BoxLib::operator>>(std::istream& is, const expect& exp)
{
    int len = exp.istr.size();
    int n = 0;
    while ( n < len )
    {
	char c;
	is >> c;
	if ( !is ) break;
	if ( c != exp.istr[n++] )
	{
	    is.putback(c);
	    break;
	}
    }
    if ( n != len )
    {
	is.clear(std::ios::badbit|is.rdstate());
	std::string msg = "expect fails to find \"" + exp.the_string() + "\"";
	BoxLib::Error(msg.c_str());
    }
    return is;
}

BoxLib::expect::expect(const char* istr_)
    : istr(istr_)
{
}

BoxLib::expect::expect(const std::string& str_)
    : istr(str_)
{
}

BoxLib::expect::expect(char c)
{
    istr += c;
}

const std::string&
BoxLib::expect::the_string() const
{
    return istr;
}
