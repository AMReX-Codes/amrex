#ifndef BL_UTILITY_H
#define BL_UTILITY_H
//
// $Id: Utility.H,v 1.51 2001-10-15 18:17:51 car Exp $
//
#include <winstd.H>

#include <iostream>
#include <string>

#ifndef WIN32
#include <sys/types.h>
#include <sys/wait.h>
#endif

#include <BLassert.H>
#include <REAL.H>

//
//@Man:
//@Memo: Useful C++ Utility Functions
/*@Doc:

  This data-less class is a poor-person's namespace of utility functions.
  Since we can't assume the existence of namespaces, and we don't
  like global functions, we put them into a data-less class as
  static member functions.
*/

namespace BoxLib
{
    /*@ManDoc: Returns the total user and system time used by
      the calling process up to the point of the call.
      If t != 0, it also stores the result in t.
    */
    double second (double* t = 0);

    /*@ManDoc: Returns the current time as the number of microseconds
               relative to program startup.  If t != 0, it also
               stores the result in t.
    */
    double wsecond (double* t = 0);
    //
    //@ManDoc: Reset start of Wall Clock Time for wsecond() to NOW.
    //
    void ResetWallClockTime ();
    //
    //@ManDoc: Return true if argument is a non-zero length string of digits.
    //
    bool is_integer (const char* str);
    //
    //@ManDoc: Returns rootNNNN where NNNN == num.
    //
    std::string Concatenate (const std::string& root,
                             int                num);

    /*@ManDoc:
      Creates the specified directories.  `path' may be either a full pathname
      or a relative pathname.  It will create all the directories in the
      pathname, if they don't already exist, so that on successful return the
      pathname refers to an existing directory.  Returns true or false
      depending upon whether or not all it was successful.  Also returns
      if `path' is NULL or "/".  `mode' is the mode passed to mkdir() for
      any directories that must be created.
    
      For example, if it is passed the string "/a/b/c/d/e/f/g", it
      will return successfully when all the directories in the pathname
      exist; i.e. when the full pathname is a valid directory.

      In a Windows environment, the path separator is a '\', so that if using
      the example given above you must pass the string  
      "\\a\\b\\c\\d\\e\\f\\g"  (Note that you must escape the backslash in a
      character string), 
    */
    bool UtilCreateDirectory (const std::string& path,
                              mode_t             mode);
    //
    //@ManDoc: Output a message and abort when couldn't create the directory.
    //
    void CreateDirectoryFailed (const std::string& dir);
    //
    //@ManDoc: Output a message and abort when couldn't open the file.
    //
    void FileOpenFailed (const std::string& dir);
    //
    //@ManDoc: Attempt to unlink the file.  Ignore any errors.
    //
    void UnlinkFile (const std::string& file);

    /*@ManDoc: Aborts after printing message indicating out-of-memory;
               i.e. operator new has failed. This is the "supported"
               set\_new\_handler() function for BoxLib applications.
     */
    void OutOfMemory ();

    /*@ManDoc:
      Mersenne Twister pseudo-random number generator.

      Generates one pseudorandom real number (double) which is
      uniformly distributed on [0,1]-interval for each call.

      Accepts any 32-bit integer as a seed -- uses 4357 as the default.

      Has a period of 2**19937.

      Mersenne Twister Home Page: http://www.math.keio.ac.jp/matumoto/emt.html

      There is also an entry point for Fortran callable as:
    
      REAL\_T rn
      call blutilrand(rn)
    */
    double Random ();

    /*@ManDoc: Set the seed of the random number generator.
    */
    void InitRandom (unsigned long seed);

    /*@ManDoc: Try to execute a shell command in a subshell.
               User is responsible for calling waitpid() on the pid
	       returned from this function.
    */
    pid_t Execute (const char* cmd);
    //
    // Helper class for Times, used in Threads, Profiling, and Timers.
    //
    class Time
    {
    public:
	Time();
	Time(long s, long n = 0);
	Time(double d);
	double as_double() const;
	long as_long() const;
	Time& operator+=(const Time&);
	Time operator+(const Time&) const;
	static Time get_time();
    private:
	long tv_sec;
	long tv_nsec;
	void normalize();
    };
    //
    // Implements a simple in place Timer:
    //
    template <double (*FCN)(double*)> class base_Timer;
    template <double (*FCN)(double*)> std::ostream& operator<<(std::ostream&, const base_Timer<FCN>&);

    typedef base_Timer<BoxLib::wsecond> WallTimer;
    typedef base_Timer<BoxLib::second>  CPUTimer;

    template <double (*FCN)(double*)>
    class base_Timer
    {
    public:
	class bad_timer;
	base_Timer();
	~base_Timer();
	void start();
	void stop();
	void reset();
	double time() const;
	double accum_time() const;
	int count() const;
	bool is_running() const;
	static double tick();
    private:
	bool running;
	double val;
	double held;
	double accum_val;
	int cnt;
    };
    //
    // The Mersenne twistor :
    //
    class mt19937
    {
    public:
	typedef unsigned long seed_type;
	explicit mt19937(seed_type seed = 4357UL);
	void rewind();
	
	double d1_value();	// [0,1] random numbers
	double d_value();	// [0,1) random numbers
	long l_value();		// [0,2^31-1] random numbers
	unsigned long u_value(); // [0,2^32-1] random numbers
    private:
	void sgenrand(unsigned long seed);
	unsigned long igenrand();
	void reload();
    private:
	enum { N = 624 };
	unsigned long init_seed;
	unsigned long mt[N];	// the array for the state vector
	int mti;		// mti==N+1 means mt[N] is not initialized
    };

    class expect;
    std::istream& operator>>(std::istream&, const expect& exp);
    
    class expect
    {
	friend std::istream& operator>>(std::istream&, const expect& exp);
    public:
	explicit expect(const std::string& str_);
	explicit expect(const char* istr_);
	explicit expect(char c);
	const std::string& the_string() const;
    private:
	std::string istr;
    };
}

//
// I'm going to document right here all the BL macros that aren't documented
// anywhere else.  Note that all these #ifdef ... #endif blocks are necessary
// to get doc++ to properly document the macros.
//

#ifdef BL_ARCH_CRAY
#undef BL_ARCH_CRAY
/*@ManDoc:
  The macro BL\_ARCH\_CRAY is defined only when compiling on a Cray
  architecture.  The make subsystem automatically sets this if it
  detects that a compilation is taking place on a Cray system.
*/
#define BL_ARCH_CRAY 1
#endif /*BL_ARCH_CRAY*/

#ifdef BL_LANG_FORT
#undef BL_LANG_FORT
/*@ManDoc:
  The macro BL\_LANG\_FORT indicates that Fortran code is being compiled.
*/
#define BL_LANG_FORT 1
#endif /*BL_LANG_FORT*/

#ifdef BL_LANG_CC
#undef BL_LANG_CC
/*@ManDoc:
  The macro BL\_LANG\_CC indicates that C++ code is being compiled.
*/
#define BL_LANG_CC 1
#endif /*BL_LANG_CC*/

#ifdef BL_FORT_USE_UNDERSCORE
#undef BL_FORT_USE_UNDERSCORE
/*@ManDoc:
  The macro BL\_FORT\_USE\_UNDERSCORE indicates that C++ code should call
  Fortran routines by appending an underscore to the name of the Fortran
  routine.  This is set automatically by the make subsystem.  It is the
  default for all BoxLib-style programs except on Cray architectures.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyz\_.
*/
#define BL_FORT_USE_UNDERSCORE 1
#endif /*BL_FORT_USE_UNDERSCORE*/

#ifdef BL_FORT_USE_UPPERCASE
#undef BL_FORT_USE_UPPERCASE
/*@ManDoc:
  The macro BL\_FORT\_USE\_UPPERCASE indicates that C++ code should call
  Fortran routines using uppercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.  Currently it
  is only set when compiling on Cray architectures.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as ABCXYZ.
*/
#define BL_FORT_USE_UPPERCASE 1
#endif /*BL_FORT_USE_UPPERCASE*/

#ifdef BL_FORT_USE_LOWERCASE
#undef BL_FORT_USE_LOWERCASE
/*@ManDoc:
  The macro BL\_FORT\_USE\_LOWERCASE indicates that C++ code should call
  Fortran routines using lowercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.  

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyx.
*/
#define BL_FORT_USE_LOWERCASE 1
#endif /*BL_FORT_USE_LOWERCASE*/

/*@ManDoc:
  BL\_IGNORE\_MAX is a macro that expands to the literal value 100000.  It is
  defined when compiling either Fortran or C++ code; i.e. when either
  BL\_LANG\_CC or BL\_LANG\_FORT is defined.  It is used in calls to
  istream::ignore() in the BoxLib code when reading in characters from an
  istream.  We use this macro instead of the more proper INT\_MAX from
  <limits.h> since at least one compiler didn't work properly when
  istream::ignore() was passed INT\_MAX.
*/
#define BL_IGNORE_MAX 100000

// Implementation of the Timer

template <double (*FCN)(double*)>
BoxLib::base_Timer<FCN>::base_Timer()
    : running(false), val(0.0), accum_val(0.0), cnt(0)
{
}

template <double (*FCN)(double*)>
BoxLib::base_Timer<FCN>::~base_Timer()
{
    BL_ASSERT( !running );
}

template <double (*FCN)(double*)>
bool
BoxLib::base_Timer<FCN>::is_running() const
{
    return running;
}

template <double (*FCN)(double*)>
inline
void
BoxLib::base_Timer<FCN>::start()
{
    BL_ASSERT( !running );
    held = FCN(0);
    running = true;
}

template <double (*FCN)(double*)>
inline
void
BoxLib::base_Timer<FCN>::stop()
{
    BL_ASSERT( running );
    val = (FCN(0) - held);
#ifndef NDEBUG
    if ( val < 0 )
    {
	std::cout << "Got a negative time " << val << "!" << std::endl;
    }
#endif
    accum_val += val;
    cnt += 1;
    running = false;
}

template <double (*FCN)(double*)>
void
BoxLib::base_Timer<FCN>::reset()
{
    BL_ASSERT( !running );
    accum_val = 0;
    cnt = 0;
}

template <double (*FCN)(double*)>
double
BoxLib::base_Timer<FCN>::accum_time() const
{
    BL_ASSERT( !running );
    return accum_val;
}

template <double (*FCN)(double*)>
double
BoxLib::base_Timer<FCN>::time() const
{
    BL_ASSERT( !running );
    return val;
}

template <double (*FCN)(double*)>
int
BoxLib::base_Timer<FCN>::count() const
{
    BL_ASSERT( !running );
    return cnt;
}

// I stole this code from someplace, but I don't know where.
template <double (*FCN)(double*)>
double
BoxLib::base_Timer<FCN>::tick()
{
    const int M = 100;
    double timesfound[M];

    //  Collect a sequence of M unique time values from the system.

    for ( int i = 0; i < M; ++i )
    {
	double t2;
	double t1 = FCN(0);
	while( ((t2 = FCN(0)) - t1) == 0 )
	{
	}
	timesfound[i] = t2;
    }

    double minDelta = timesfound[1] - timesfound[0];
    for ( int i = 2; i < M; i++ )
    {
	minDelta = std::min(minDelta, std::max(timesfound[i]-timesfound[i-1], 0.0));
    }
    return minDelta;
}

template <double (*FCN)(double*)>
std::ostream& operator<<(std::ostream& os, const BoxLib::base_Timer<FCN>& bt)
{
    return os << "["
	      << bt.accum_time() << "/"
	      << bt.count()
	      << "]";
}

#endif /*BL_UTILITY_H*/
