//
// $Id: Utility.cpp,v 1.84 2009-12-21 20:36:20 lijewski Exp $
//

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <BLFort.H>
#include <REAL.H>
#include <BoxLib.H>
#include <Utility.H>
#include <BLassert.H>
#include <Profiler.H>

#if defined(BL_XT3) && !defined(BL_USE_MPI)
#include <unistd.h>
#endif

#ifdef BL_BGL
#include <ParallelDescriptor.H>
#endif

#ifdef WIN32
#include <direct.h>
#define mkdir(a,b) _mkdir((a))
const char* path_sep_str = "\\";
#else
const char* path_sep_str = "/";
#endif

#if !defined(WIN32) && !defined(BL_XT3)

#include <sys/types.h>
#include <sys/times.h>
#ifdef BL_AIX
#undef _XOPEN_SOURCE_EXTENDED
#define _XOPEN_SOURCE_EXTENDED 1
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

#elif defined(WIN32)

// minimum requirement of WindowsNT
#include <windows.h>

namespace
{
double rate;
bool inited = false;
LONGLONG
get_initial_wall_clock_time()
{
    LARGE_INTEGER li;
    QueryPerformanceFrequency(&li);
    rate = 1.0/li.QuadPart;
    QueryPerformanceCounter(&li);
    inited = true;
    return li.QuadPart;
}
LONGLONG BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();
}
double
BoxLib::wsecond(double* rslt)
{
    BL_ASSERT( inited );
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    double result = double(li.QuadPart-BL_Initial_Wall_Clock_Time)*rate;
    if ( rslt ) *rslt = result;
    return result;
}

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

#elif defined(BL_XT3) && defined(BL_USE_MPI)

#include <catamount/dclock.h>
#include <unistd.h>

static
double
get_initial_wall_clock_time ()
{
    return dclock();
}

//
// Attempt to guarantee wsecond() gets initialized on program startup.
//
double BL_Initial_Wall_Clock_Time = get_initial_wall_clock_time();

double
BoxLib::wsecond (double* t_)
{
    double t = dclock() - BL_Initial_Wall_Clock_Time;
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

#endif

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

std::string
BoxLib::Concatenate (const std::string& root,
                     int                num,
                     int                mindigits)
{
    std::string result = root;
    char buf[32];
    sprintf(buf, "%0*d",  mindigits, num);
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
BoxLib::UtilCreateDirectory (const std::string& path,
                             int)
#else
BoxLib::UtilCreateDirectory (const std::string& path,
                             mode_t         mode)
#endif
{
    BL_PROFILE("BoxLib::UtilCreateDirectory()");
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
BoxLib::CreateDirectoryFailed (const std::string& dir)
{
    std::string msg("Couldn't create directory: ");
    msg += dir;
    BoxLib::Error(msg.c_str());
}

void
BoxLib::FileOpenFailed (const std::string& file)
{
    std::string msg("Couldn't open file: ");
    msg += file;
    BoxLib::Error(msg.c_str());
}

void
BoxLib::UnlinkFile (const std::string& file)
{
    unlink(file.c_str());
}

void
BoxLib::OutOfMemory ()
{
#ifdef BL_BGL
    ParallelDescriptor::Abort(12);
#else
    BoxLib::Error("Sorry, out of memory, bye ...");
#endif
}

#if 0
#ifdef WIN32
pid_t
BoxLib::Execute (const char* cmd)
{
  BoxLib::Error("Execute failed!");
  return -1;
}
#else
//extern "C" pid_t fork();

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

//
// initializing with a NONZERO seed.
//
void
BoxLib::mt19937::sgenrand(unsigned long seed)
{
    mt[0]= seed & 0xffffffffUL;
    for ( mti=1; mti<N; mti++ ) 
    {
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30L)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;       /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void 
BoxLib::mt19937::sgenrand(unsigned long init_key[], int key_length)
{
    int i, j, k;
    sgenrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for ( ; k; k-- ) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for ( k=N-1; k; k-- ) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

void
BoxLib::mt19937::reload()
{
    unsigned long y;
    int kk;

    const int M = 397;

#define MATRIX_A    0x9908B0DFUL // Constant vector a
#define UPPER_MASK  0x80000000UL // Most significant w-r bits
#define LOWER_MASK  0x7FFFFFFFUL // least significant r bits
    //
    // mag01[x] = x * MATRIX_A  for x=0,1
    //
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    for ( kk=0; kk<N-M; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+M] ^ (y >> 1L) ^ mag01[y & 0x1UL];
    }
    for ( ; kk<N-1; kk++ )
    {
	y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	mt[kk] = mt[kk+(M-N)] ^ (y >> 1L) ^ mag01[y & 0x1UL];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1L) ^ mag01[y & 0x1UL];

    mti = 0;

#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
}

unsigned long
BoxLib::mt19937::igenrand()
{
    //
    // Generate N words at one time.
    //
    if ( mti >= N ) reload();

    unsigned long y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

BoxLib::mt19937::mt19937(unsigned long seed)
    :
    init_seed(seed),
    mti(N)
{
    sgenrand(seed);
}

BoxLib::mt19937::mt19937 (unsigned long seed_array[], int len)
{
    sgenrand(seed_array, len);
}

void
BoxLib::mt19937::rewind()
{
    sgenrand(init_seed);
}

//
// [0,1] random numbers
//
double
BoxLib::mt19937::d_value()
{
    return double(igenrand()) * (1.0/4294967295.0);  // divided by 2^32-1
}

//
// [0,1) random numbers
//
double
BoxLib::mt19937::d1_value()
{
    return double(igenrand()) * (1.0/4294967296.0);  // divided by 2^32
}

//
// (0,1) random numbers
//
double
BoxLib::mt19937::d2_value()
{
    return (double(igenrand()) + .5) * (1.0/4294967296.0);  // divided by 2^32
}

long
BoxLib::mt19937::l_value()
{
    return (long)(igenrand()>>1);
}

unsigned long
BoxLib::mt19937::u_value()
{
    return igenrand();
}

void
BoxLib::mt19937::save (Array<unsigned long>& state) const
{
    state.resize(N+2);
    state[0] = init_seed;
    for (int i = 0; i < N; i++)
        state[i+1] = mt[i];
    state[N+1] = mti;
}

void
BoxLib::mt19937::restore (const Array<unsigned long>& state)
{
    if (state.size() != N+2)
        BoxLib::Error("mt19937::restore(): incorrectly sized state vector");

    init_seed = state[0];
    for (int i = 0; i < N; i++)
        mt[i] = state[i+1];
    mti = state[N+1];

    if (mti < 0 || mti > N)
        BoxLib::Error("mt19937::restore(): mti out-of-bounds");
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

void
BoxLib::SaveRandomState (Array<unsigned long>& state)
{
    the_generator.save(state);
}

void
BoxLib::RestoreRandomState (const Array<unsigned long>& state)
{
    the_generator.restore(state);
}

//
// Fortran entry points for BoxLib::Random().
//

BL_FORT_PROC_DECL(BLUTILINITRAND,blutilinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    BoxLib::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLINITRAND,blinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    BoxLib::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLUTILRAND,blutilrand)(Real* rn)
{
    *rn = BoxLib::Random();
}

//
// Lower tail quantile for standard normal distribution function.
//
// This function returns an approximation of the inverse cumulative
// standard normal distribution function.  I.e., given P, it returns
// an approximation to the X satisfying P = Pr{Z <= X} where Z is a
// random variable from the standard normal distribution.
//
// The algorithm uses a minimax approximation by rational functions
// and the result has a relative error whose absolute value is less
// than 1.15e-9.
//
// Author:      Peter J. Acklam
// Time-stamp:  2002-06-09 18:45:44 +0200
// E-mail:      jacklam@math.uio.no
// WWW URL:     http://www.math.uio.no/~jacklam
//
// C implementation adapted from Peter's Perl version
//

double
BoxLib::InvNormDist (double p)
{
    if (p <= 0 || p >= 1)
        BoxLib::Error("BoxLib::InvNormDist(): p MUST be in (0,1)");
    //
    // Coefficients in rational approximations.
    //
    static const double a[6] =
    {
	-3.969683028665376e+01,
        2.209460984245205e+02,
	-2.759285104469687e+02,
        1.383577518672690e+02,
	-3.066479806614716e+01,
        2.506628277459239e+00
    };
    static const double b[5] =
    {
	-5.447609879822406e+01,
        1.615858368580409e+02,
	-1.556989798598866e+02,
        6.680131188771972e+01,
	-1.328068155288572e+01
    };
    static const double c[6] =
    {
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00
    };
    static const double d[4] =
    {
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
    };

    static const double lo = 0.02425;
    static const double hi = 0.97575;

    double x;

    if (p < lo)
    {
        //
        // Rational approximation for lower region.
        //
        double q = std::sqrt(-2*std::log(p));

        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else if (p > hi)
    {
        //
        // Rational approximation for upper region.
        //
        double q = std::sqrt(-2*std::log(1-p));

        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else
    {
        //
        // Rational approximation for central region.
        //
        double q = p - 0.5;
        double r = q*q;

        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }

    return x;
}

BL_FORT_PROC_DECL(BLINVNORMDIST,blinvnormdist)(Real* result)
{
    //
    // Get a random number in (0,1);
    //
    double val = the_generator.d2_value();

    *result = BoxLib::InvNormDist(val);
}

//
//****************************************************************************80
//
//  Purpose:
//
//    InvNormDistBest inverts the standard normal CDF.
//
//  Discussion:
//
//    The result is accurate to about 1 part in 10**16.
//
//  Modified:
//
//    27 December 2004
//
//  Author:
//
//    Original FORTRAN77 version by Michael Wichura.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Wichura,
//    The Percentage Points of the Normal Distribution,
//    Algorithm AS 241,
//    Applied Statistics,
//    Volume 37, Number 3, pages 477-484, 1988.
//
//  Parameters:
//
//    Input, double P, the value of the cumulative probability 
//    densitity function.  0 < P < 1.  Fails if P is outside this range.
//
//    Output, the normal deviate value with the property that the
//    probability of a standard normal deviate being less than or equal
//    to this value is P.
//

double
BoxLib::InvNormDistBest (double p)

{
  static const double a[8] =
  { 
      3.3871328727963666080,     1.3314166789178437745e+2,
      1.9715909503065514427e+3,  1.3731693765509461125e+4,
      4.5921953931549871457e+4,  6.7265770927008700853e+4,
      3.3430575583588128105e+4,  2.5090809287301226727e+3
  };
  static const double b[8] =
  {
      1.0,                       4.2313330701600911252e+1,
      6.8718700749205790830e+2,  5.3941960214247511077e+3,
      2.1213794301586595867e+4,  3.9307895800092710610e+4,
      2.8729085735721942674e+4,  5.2264952788528545610e+3
  };
  static const double c[8] =
  {
      1.42343711074968357734,     4.63033784615654529590,
      5.76949722146069140550,     3.64784832476320460504,
      1.27045825245236838258,     2.41780725177450611770e-1,
      2.27238449892691845833e-2,  7.74545014278341407640e-4
  };

  static const double d[8] =
  {
      1.0,                        2.05319162663775882187,
      1.67638483018380384940,     6.89767334985100004550e-1,
      1.48103976427480074590e-1,  1.51986665636164571966e-2,
      5.47593808499534494600e-4,  1.05075007164441684324e-9
  };
  static const double e[8] =
  {
      6.65790464350110377720,     5.46378491116411436990,
      1.78482653991729133580,     2.96560571828504891230e-1,
      2.65321895265761230930e-2,  1.24266094738807843860e-3,
      2.71155556874348757815e-5,  2.01033439929228813265e-7
  };
  static const double f[8] =
  {
      1.0,                        5.99832206555887937690e-1,
      1.36929880922735805310e-1,  1.48753612908506148525e-2,
      7.86869131145613259100e-4,  1.84631831751005468180e-5,
      1.42151175831644588870e-7,  2.04426310338993978564e-15
  };

  static const double CONST1 = 0.180625;
  static const double CONST2 = 1.6;
  static const double SPLIT1 = 0.425;
  static const double SPLIT2 = 5.0;

  double r, value;

  if (p <= 0 || p >= 1)
      BoxLib::Error("InvNormDistBest(): p MUST be in (0,1)");

  double q = p - 0.5;

  if ( std::fabs ( q ) <= SPLIT1 )
  {
      r = CONST1 - q * q;

      double num = 0.0, den = 0.0;

      for (int i = 7; 0 <= i; i-- )
      {
          num = num * r + a[i];
          den = den * r + b[i];
      }

      value = q * num / den;
  }
  else
  {
      r = ( ( q < 0.0 ) ? p : (1.0 - p) );

      r = std::sqrt ( -std::log ( r ) );

      if ( r <= SPLIT2 )
      {
          r = r - CONST2;

          double num = 0.0, den = 0.0;

          for (int i = 7; 0 <= i; i-- )
          {
              num = num * r + c[i];
              den = den * r + d[i];
          }

          value = num / den;
      }
      else
      {
          r = r - SPLIT2;

          double num = 0.0, den = 0.0;

          for (int i = 7; 0 <= i; i-- )
          {
              num = num * r + e[i];
              den = den * r + f[i];
          }

          value = num / den;
      }

      if ( q < 0.0 ) value = -value;
  }

  return value;
}

BL_FORT_PROC_DECL(BLINVNORMDISTBEST,blinvnormdistbest)(Real* result)
{
    //
    // Get a random number in (0,1);
    //
    double val = the_generator.d2_value();

    *result = BoxLib::InvNormDistBest(val);
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
