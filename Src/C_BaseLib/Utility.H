
#ifndef BL_UTILITY_H
#define BL_UTILITY_H

#include <winstd.H>

#include <iostream>
#include <string>
#include <map>

#ifndef WIN32
#include <sys/types.h>
#include <unistd.h>
#if !defined(COLUMBIA)
#include <sys/wait.h>
#endif
#endif

#include <BLassert.H>
#include <REAL.H>
#include <Array.H>
#include <Box.H>
#include <BoxArray.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>

//
// Useful C++ Utility Functions
//
// This data-less class is a poor-person's namespace of utility functions.
// Since we can't assume the existence of namespaces, and we don't
// like global functions, we put them into a data-less class as
// static member functions.
//

namespace BoxLib
{
    //
    // Returns the total user and system time used by
    // the calling process up to the point of the call.
    // If t != 0, it also stores the result in t.
    //
    double second (double* t = 0);
    //
    // Returns the current time as the number of microseconds
    // relative to program startup.  If t != 0, it also
    // stores the result in t.
    //
    double wsecond (double* t = 0);
    //
    // Reset start of Wall Clock Time for wsecond() to NOW.
    //
    void ResetWallClockTime ();
    //
    // Return true if argument is a non-zero length string of digits.
    //
    bool is_integer (const char* str);
    //
    // Splits "instr" into separate pieces based on "separators".
    //
    const std::vector<std::string>& Tokenize (const std::string& instr,
                                              const std::string& separators);
    //
    // Returns rootNNNN where NNNN == num.
    //
    std::string Concatenate (const std::string& root,
                             int                num,
                             int                mindigits = 5);
    //
    //  Creates the specified directories.  path may be either a full pathname
    //  or a relative pathname.  It will create all the directories in the
    //  pathname, if they don't already exist, so that on successful return the
    //  pathname refers to an existing directory.  Returns true or false
    //  depending upon whether or not it was successful.  Also returns true
    //  if path is NULL or "/".  mode is the mode passed to mkdir() for
    //  any directories that must be created (for example:  0755).
    //  verbose will print out the directory creation steps.
    // 
    //  For example, if it is passed the string "/a/b/c/d/e/f/g", it
    //  will return successfully when all the directories in the pathname
    //  exist; i.e. when the full pathname is a valid directory.
    //
    //  In a Windows environment, the path separator is a '\', so that if using
    //  the example given above you must pass the string  
    //  "\\a\\b\\c\\d\\e\\f\\g"  (Note that you must escape the backslash in a
    //  character string), 
    //
    //  Only the last mkdir return value is checked for success as errno may not
    //  be set to EEXIST if a directory exists but mkdir has other reasons
    //  to fail such as part of the path being a read-only filesystem (EROFS).
    //  If this function fails, it will print out an error stack.
    //
    bool UtilCreateDirectory (const std::string& path,
                              mode_t             mode,
			      bool               verbose = false);
    //
    // Output a message and abort when couldn't create the directory.
    //
    void CreateDirectoryFailed (const std::string& dir);
    //
    // Output a message and abort when couldn't open the file.
    //
    void FileOpenFailed (const std::string& dir);
    //
    // Attempt to unlink the file.  Ignore any errors.
    //
    void UnlinkFile (const std::string& file);
    //
    // Check if a file already exists.
    //   Return true if the filename is an existing file, directory,
    //   or link.  For links, this operates on the link and not what
    //   the link points to.
    //
    bool FileExists(const std::string &filename);
    //
    // Create a (probably) unique string
    //
    std::string UniqueString();
    //
    // Create a new directory, renaming the old one if it exists
    //
    void UtilCreateCleanDirectory (const std::string &path,
                                   bool callbarrier = true);
    //
    // Rename a current directory if it exists
    //
    void UtilRenameDirectoryToOld (const std::string &path,
                                   bool callbarrier = true);
    //
    // Aborts after printing message indicating out-of-memory;
    // i.e. operator new has failed. This is the "supported"
    // set_new_handler() function for BoxLib applications.
    //
    void OutOfMemory ();
    /*
      This function returns an approximation of the inverse cumulative
      standard normal distribution function.  I.e., given P, it returns
      an approximation to the X satisfying P = Pr{Z <= X} where Z is a
      random variable from the standard normal distribution.

      The algorithm uses a minimax approximation by rational functions
      and the result has a relative error whose absolute value is less
      than 1.15e-9.

      Author:      Peter J. Acklam
      Time-stamp:  2002-06-09 18:45:44 +0200
      E-mail:      jacklam@math.uio.no
      WWW URL:     http://www.math.uio.no/~jacklam

      "p" MUST be in the open interval (0,1).

      Entry points for Fortran:

        REAL_T val
        call blinvnormdist(val)

      Internally, the Fortran entry point calls a static Mersenne Twister
      object (the same one called by blutilrand()) to get a random number in
      the open interval (0,1), and then sets "val" to the result of calling
      BoxLib::InvNormDist() with that random number.
    */
    double InvNormDist (double p);
    /*
      This function returns an approximation of the inverse cumulative
      standard normal distribution function.  I.e., given P, it returns
      an approximation to the X satisfying P = Pr{Z <= X} where Z is a
      random variable from the standard normal distribution.

      Original FORTRAN77 version by Michael Wichura.

      Michael Wichura,
      The Percentage Points of the Normal Distribution,
      Algorithm AS 241,
      Applied Statistics,
      Volume 37, Number 3, pages 477-484, 1988.

      Our version is based on the C++ version by John Burkardt.

      The algorithm uses a minimax approximation by rational functions
      and the result is good to roughly machine precision.  This routine
      is roughly 30% more costly than InvNormDist() above.

      "p" MUST be in the open interval (0,1).

      Entry points for Fortran:
    
        REAL_T val
        call blinvnormdistbest(val)

      Internally, the Fortran entry point calls a static Mersenne Twister
      object (the same one called by blutilrand()) to get a random number in
      the open interval (0,1), and then sets "val" to the result of calling
      BoxLib::InvNormDist() with that random number.
    */
    double InvNormDistBest (double p);
    /*
      Mersenne Twister pseudo-random number generator.

      Generates one pseudorandom real number (double) which is
      uniformly distributed on [0,1]-interval for each call.

      Accepts any 32-bit integer as a seed -- uses 4357 as the default.

      Has a period of 2**19937.

      Mersenne Twister Home Page: http://www.math.keio.ac.jp/matumoto/emt.html

      There is also an entry point for Fortran callable as:
    
      REAL_T rn
      call blutilrand(rn)

      Internally, blutilrand() calls a static Mersenne Twister oject (the
      same one used by BoxLib::Random()) to get a value in [0,1] and then
      sets "rn" to that value.
    */
    double Random (); // [0,1]
    double Random1 (); // [0,1)
    double Random2 (); // (0,1)
    unsigned long Random_int(unsigned long n); // [0,n-1], where n<=2^32-1
    /* Set the seed of the random number generator.

      There is also an entry point for Fortran callable as:

      INTEGER seed
      call blutilinitrand(seed)

      or

      INTEGER seed
      call blinitrand(seed)
    */
    void InitRandom (unsigned long seed);
    void InitRandom (unsigned long seed, int numprocs);

    void ResetRandomSeed(unsigned long seed);
    //
    // Save and restore random state.
    //
    // state.size() == 626 on return from Save & on entry to Restore.
    //
    void SaveRandomState (Array<unsigned long>& state);

    int sizeofRandomState ();

    void RestoreRandomState (const Array<unsigned long>& state);
    //
    // Create a unique subset of random numbers from a pool
    //   of integers in the range [0, poolSize - 1]
    //   the set will be in the order they are found
    //   setSize must be <= poolSize
    //   uSet will be resized to setSize
    //   if you want all processors to have the same set,
    //   call this on one processor and broadcast the array
    //
    void UniqueRandomSubset (Array<int> &uSet, int setSize, int poolSize,
                             bool printSet = false);

    //
    // Fill binCounts with the number of items for each bin.
    //   The remainder items are randomly assigned.
    //
    void NItemsPerBin (int totalItems, Array<int> &binCounts);
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

	explicit mt19937 (seed_type seed = 4357UL);
        mt19937 (seed_type seed, int numprocs);
	mt19937 (seed_type array[], int array_len);
	void rewind();
	void reset(unsigned long seed);

	double d_value  (); // [0,1] random numbers
	double d1_value (); // [0,1) random numbers
	double d2_value (); // (0,1) random numbers

	long          l_value (); // [0,2^31-1] random numbers
	unsigned long u_value (); // [0,2^32-1] random numbers

        void save (Array<unsigned long>& state) const;
        int  RNGstatesize() const;
        void restore (const Array<unsigned long>& state);
    private:
	enum { N = 624 };
        static unsigned long init_seed;
        static unsigned long mt[N];	// the array for the state vector
        static int mti;		// mti==N+1 means mt[N] is not initialized
#ifdef _OPENMP
#pragma omp threadprivate(init_seed,mt,mti)
#endif
    private:
	void sgenrand (unsigned long seed);
	void sgenrand (seed_type seed_array[], int len);
	unsigned long igenrand ();
	void reload ();
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

    class StreamRetry
    {
      public:
        StreamRetry(std::ostream &os, const std::string &suffix,
	            const int maxtries);
        StreamRetry(const std::string &filename,
	            const bool abortonretryfailure,
	            const int maxtries);
	bool TryOutput();
	bool TryFileOutput();
        static int NStreamErrors()       { return nStreamErrors; }
        static void ClearStreamErrors()  { nStreamErrors = 0;    }

      private:
        int tries, maxTries;
	bool abortOnRetryFailure;
	std::string fileName;
	std::ostream &sros;
	std::ostream::pos_type spos;
	std::string suffix;

	static int nStreamErrors;
    };

    Array<char> SerializeStringArray(const Array<std::string> &stringArray);
    Array<std::string> UnSerializeStringArray(const Array<char> &charArray);
    void SyncStrings(const Array<std::string> &localStrings,
                     Array<std::string> &syncedStrings, bool &alreadySynced);

    //
    // Memory Usage Counting
    //
    // For gcc, there are extra 32 bytes for each map node.
    // For others, this number may be different.
    static const long gcc_map_node_extra_bytes = 32L;
    template <typename T> long bytesOf (const std::vector<T>& v);
    template <typename Key, typename T, class Compare> long bytesOf (const std::map<Key,T,Compare>& m);

    template<class T> void BroadcastArray(Array<T> &aT, int myLocalId, int rootId, const MPI_Comm &localComm);
    void BroadcastBox(Box &bB, int myLocalId, int rootId, const MPI_Comm &localComm);
    void BroadcastBoxArray(BoxArray &bBA, int myLocalId, int rootId, const MPI_Comm &localComm);
    void BroadcastDistributionMapping(DistributionMapping &dM, int sentinelProc,
                                      int myLocalId, int rootId, const MPI_Comm &localComm,
				      bool addToCache);
    const unsigned int msps(1000000);
    void USleep(double sleepsec);

}

template<class T> void BoxLib::BroadcastArray(Array<T> &aT, int myLocalId, int rootId, const MPI_Comm &localComm)
{
  int aT_Size(-2);
  if(myLocalId == rootId) {
    aT_Size = aT.size();
  }
  ParallelDescriptor::Bcast(&aT_Size, 1, rootId, localComm);
  BL_ASSERT(aT_Size >= 0);
  if(myLocalId != rootId) {
    aT.resize(aT_Size);
  }
  if(aT_Size > 0) {
    ParallelDescriptor::Bcast(aT.dataPtr(), aT.size(), rootId, localComm);
  }
}


//
// I'm going to document right here all the BL macros that aren't documented
// anywhere else.  Note that all these #ifdef ... #endif blocks are necessary
// to get doc++ to properly document the macros.
//

#ifdef BL_LANG_FORT
#undef BL_LANG_FORT
/*
  The macro BL_LANG_FORT indicates that Fortran code is being compiled.
*/
#define BL_LANG_FORT 1
#endif /*BL_LANG_FORT*/

#ifdef BL_FORT_USE_UNDERSCORE
#undef BL_FORT_USE_UNDERSCORE
/*
  The macro BL_FORT_USE_UNDERSCORE indicates that C++ code should call
  Fortran routines by appending an underscore to the name of the Fortran
  routine.  This is set automatically by the make subsystem.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyz_.
*/
#define BL_FORT_USE_UNDERSCORE 1
#endif /*BL_FORT_USE_UNDERSCORE*/

#ifdef BL_FORT_USE_UPPERCASE
#undef BL_FORT_USE_UPPERCASE
/*
  The macro BL_FORT_USE_UPPERCASE indicates that C++ code should call
  Fortran routines using uppercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as ABCXYZ.
*/
#define BL_FORT_USE_UPPERCASE 1
#endif /*BL_FORT_USE_UPPERCASE*/

#ifdef BL_FORT_USE_LOWERCASE
#undef BL_FORT_USE_LOWERCASE
/*
  The macro BL_FORT_USE_LOWERCASE indicates that C++ code should call
  Fortran routines using lowercase letters for all the letters in the
  routine.  This is set automatically by the make subsystem.  

  For example, if the Fortran routine is named abcxyx, then it will
  be called in C++ code as abcxyx.
*/
#define BL_FORT_USE_LOWERCASE 1
#endif /*BL_FORT_USE_LOWERCASE*/

/*
  BL_IGNORE_MAX is a macro that expands to the literal value 100000.  It is
  defined when compiling either Fortran or C++ code; i.e. when either
  BL_LANG_CC or BL_LANG_FORT is defined.  It is used in calls to
  istream::ignore() in the BoxLib code when reading in characters from an
  istream.  We use this macro instead of the more proper INT_MAX from
  <limits.h> since at least one compiler didn't work properly when
  istream::ignore() was passed INT_MAX.
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

template <typename T> 
long 
BoxLib::bytesOf (const std::vector<T>& v)
{
    return sizeof(v) + v.capacity()*sizeof(T);
}

template <typename Key, typename T, class Compare>
long
BoxLib::bytesOf (const std::map<Key,T,Compare>& m)
{
    return sizeof(m) + m.size()*(sizeof(Key)+sizeof(T)+gcc_map_node_extra_bytes);
}

#endif /*BL_UTILITY_H*/
