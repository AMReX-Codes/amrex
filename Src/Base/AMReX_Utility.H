#ifndef BL_UTILITY_H
#define BL_UTILITY_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <type_traits>

#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParallelDescriptor.H>



namespace amrex
{
/**
* \brief Useful C++ Utility Functions
*/

    //! This shuts up the compiler about unused variables
    template<class T> void ignore_unused( const T& ) { }

    //! Return true if argument is a non-zero length string of digits.
    bool is_integer (const char* str);
    //! Splits "instr" into separate pieces based on "separators".
    const std::vector<std::string>& Tokenize (const std::string& instr,
                                              const std::string& separators);

    //! Converts all characters of the string into lower or uppercase based on std::locale 
    std::string toLower (std::string s);
    std::string toUpper (std::string s);

    //! Returns rootNNNN where NNNN == num.
    std::string Concatenate (const std::string& root,
                             int                num,
                             int                mindigits = 5);
    /**
    *  \brief Creates the specified directories.  path may be either a full pathname
    *  or a relative pathname.  It will create all the directories in the
    *  pathname, if they don't already exist, so that on successful return the
    *  pathname refers to an existing directory.  Returns true or false
    *  depending upon whether or not it was successful.  Also returns true
    *  if path is NULL or "/".  mode is the mode passed to mkdir() for
    *  any directories that must be created (for example:  0755).
    *  verbose will print out the directory creation steps.
    * 
    *  For example, if it is passed the string "/a/b/c/d/e/f/g", it
    *  will return successfully when all the directories in the pathname
    *  exist; i.e. when the full pathname is a valid directory.
    *
    *  In a Windows environment, the path separator is a '\', so that if using
    *  the example given above you must pass the string  
    *  "\\a\\b\\c\\d\\e\\f\\g"  (Note that you must escape the backslash in a
    *  character string), 
    *
    *  Only the last mkdir return value is checked for success as errno may not
    *  be set to EEXIST if a directory exists but mkdir has other reasons
    *  to fail such as part of the path being a read-only filesystem (EROFS).
    *  If this function fails, it will print out an error stack.
    */
    bool UtilCreateDirectory (const std::string& path,
                              mode_t             mode,
			      bool               verbose = false);
    //! Output a message and abort when couldn't create the directory.
    void CreateDirectoryFailed (const std::string& dir);
    //! Output a message and abort when couldn't open the file.
    void FileOpenFailed (const std::string& dir);
    //! Attempt to unlink the file.  Ignore any errors.
    void UnlinkFile (const std::string& file);
    /**
    * \brief Check if a file already exists.
    *   Return true if the filename is an existing file, directory,
    *   or link.  For links, this operates on the link and not what
    *   the link points to.
    */
    bool FileExists(const std::string &filename);
    //! Create a (probably) unique string
    std::string UniqueString();
    //! Create a new directory, renaming the old one if it exists
    void UtilCreateCleanDirectory (const std::string &path,
                                   bool callbarrier = true);

    ///
    /**
       Create a new directory, removing old one if it exists.
       This will only work on unix systems, as it has a system call.
     */
    void UtilCreateDirectoryDestructive(const std::string &path,
                                        bool callbarrier = true);

    //! Rename a current directory if it exists
    void UtilRenameDirectoryToOld (const std::string &path,
                                   bool callbarrier = true);
    /**
    * \brief Aborts after printing message indicating out-of-memory;
    * i.e. operator new has failed. This is the "supported"
    * set_new_handler() function for AMReX applications.
    */
    void OutOfMemory ();
    /**
    *  \brief This function returns an approximation of the inverse cumulative
    *  standard normal distribution function.  I.e., given P, it returns
    *  an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    *  random variable from the standard normal distribution.

    *  The algorithm uses a minimax approximation by rational functions
    *  and the result has a relative error whose absolute value is less
    *  than 1.15e-9.

    *  \author      Peter J. Acklam
    *  Time-stamp:  2002-06-09 18:45:44 +0200
    *  E-mail:      jacklam@math.uio.no
    *  WWW URL:     http://www.math.uio.no/~jacklam

    *  "p" MUST be in the open interval (0,1).

    *  Entry points for Fortran:

    *    REAL_T val
    *    call blinvnormdist(val)

    *  Internally, the Fortran entry point calls a static Mersenne Twister
    *  object (the same one called by blutilrand()) to get a random number in
    *  the open interval (0,1), and then sets "val" to the result of calling
    *  amrex::InvNormDist() with that random number.
    */
    double InvNormDist (double p);
    /**
    * \brief  This function returns an approximation of the inverse cumulative
    *  standard normal distribution function.  I.e., given P, it returns
    *  an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    *  random variable from the standard normal distribution.
    *
    *  Original FORTRAN77 version by Michael Wichura.
    *
    *  Michael Wichura,
    *  The Percentage Points of the Normal Distribution,
    *  Algorithm AS 241,
    *  Applied Statistics,
    *  Volume 37, Number 3, pages 477-484, 1988.
    *
    *  Our version is based on the C++ version by John Burkardt.

    *  The algorithm uses a minimax approximation by rational functions
    *  and the result is good to roughly machine precision.  This routine
    *  is roughly 30% more costly than InvNormDist() above.
    *
    *  "p" MUST be in the open interval (0,1).
    *
    *  Entry points for Fortran:
    * 
    *    REAL_T val
    *    call blinvnormdistbest(val)
    *
    *  Internally, the Fortran entry point calls a static Mersenne Twister
    *  object (the same one called by blutilrand()) to get a random number in
    *  the open interval (0,1), and then sets "val" to the result of calling
    *  amrex::InvNormDist() with that random number.
    */
    double InvNormDistBest (double p);

    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) from a normal
    *  distribution with mean 'mean' and standard deviation 'stddev'.
    *
    */
    AMREX_GPU_HOST_DEVICE double RandomNormal (double mean, double stddev);

    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) which is
    *  uniformly distributed on [0,1)-interval for each call.
    *
    */
    AMREX_GPU_HOST_DEVICE double Random (); // [0,1)

    /**
    * \brief Generate a psuedo-random int using C++11's mt19937.
    *
    *  Generates one pseudorandom unsigned long integer which is
    *  uniformly distributed on [0,n-1]-interval for each call.
    *
    */
    unsigned long Random_int(unsigned long n); // [0,n-1]

    /** \brief Set the seed of the random number generator.
    *
    *  There is also an entry point for Fortran callable as:
    *
    *  INTEGER seed
    *  call blutilinitrand(seed)
    *
    *  or
    *
    *  INTEGER seed
    *  call blinitrand(seed)
    */
    void InitRandom (unsigned long seed, int nprocs=ParallelDescriptor::NProcs());

    void ResetRandomSeed(unsigned long seed);
    /**
    * \brief Save and restore random state.
    *
    */
    void SaveRandomState(std::ostream& os);

    void RestoreRandomState(std::istream& is, int nthreads_old, int nstep_old);

    /**
    * \brief Create a unique subset of random numbers from a pool
    *   of integers in the range [0, poolSize - 1]
    *   the set will be in the order they are found
    *   setSize must be <= poolSize
    *   uSet will be resized to setSize
    *   if you want all processors to have the same set,
    *   call this on one processor and broadcast the array
    */
    void UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                             bool printSet = false);

    /**
    * \brief Fill binCounts with the number of items for each bin.
    *   The remainder items are randomly assigned.
    */
    void NItemsPerBin (int totalItems, Vector<int> &binCounts);

    /* \brief cumulative refinement ratio between levels */
    int CRRBetweenLevels(int fromlevel, int tolevel,
                         const Vector<int> &refratios);

    /**
    * \brief Set the seeds of the random number generator for each CUDA thread .
    *  An in-built CUDA library, curand_init is used here.
    *  The pseusorandom sequence currently implemented is obtained 
    *  from the XORWOW and MRG32k3a generator
    */
    void InitRandSeedOnDevice(int N); 
   

    /**
    * \brief Resize seed array and copy address to global symbol 
    */
    void ResizeRandomSeed(int N);

    void CheckSeedArraySizeAndResize (int N);

    void DeallocateRandomSeedDevArray();

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

    Vector<char> SerializeStringArray(const Vector<std::string> &stringArray);
    Vector<std::string> UnSerializeStringArray(const Vector<char> &charArray);
    void SyncStrings(const Vector<std::string> &localStrings,
                     Vector<std::string> &syncedStrings, bool &alreadySynced);

    //
    // Memory Usage Counting
    //
    // For gcc, there are extra 32 bytes for each map node.
    // For others, this number may be different.
    static const long gcc_map_node_extra_bytes = 32L;
    template <typename T> long bytesOf (const std::vector<T>& v);
    template <typename Key, typename T, class Compare> long bytesOf (const std::map<Key,T,Compare>& m);

    void BroadcastBool(bool &bBool, int myLocalId, int rootId, const MPI_Comm &localComm);

    void BroadcastString(std::string &bStr, int myLocalId, int rootId, const MPI_Comm &localComm);
    void BroadcastStringArray(Vector<std::string> &bSA, int myLocalId, int rootId, const MPI_Comm &localComm);

    template<class T> void BroadcastArray(Vector<T> &aT, int myLocalId, int rootId, const MPI_Comm &localComm);
    void USleep(double sleepsec);


    using MaxResSteadyClock = std::conditional<std::chrono::high_resolution_clock::is_steady,
                                               std::chrono::high_resolution_clock,
                                               std::chrono::steady_clock>::type;
    double second () noexcept;

    template<typename T> void hash_combine (uint64_t & seed, const T & val) noexcept;
    template<typename T> uint64_t hash_vector (const Vector<T> & vec, uint64_t seed = 0xDEADBEEFDEADBEEF) noexcept;

    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE const T& min (const T& a, const T& b) noexcept
    {
        return (b < a) ? b : a;
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE const T& min (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return min(min(a,b),c...);
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE const T& max (const T& a, const T& b) noexcept
    {
        return (a < b) ? b : a;
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE const T& max (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return max(max(a,b),c...);
    }
}

template<class T> void amrex::BroadcastArray(Vector<T> &aT, int myLocalId, int rootId, const MPI_Comm &localComm)
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
  istream::ignore() in the AMReX code when reading in characters from an
  istream.  We use this macro instead of the more proper INT_MAX from
  <limits> since at least one compiler didn't work properly when
  istream::ignore() was passed INT_MAX.
*/
#define BL_IGNORE_MAX 100000

template <typename T> 
long 
amrex::bytesOf (const std::vector<T>& v)
{
    return sizeof(v) + v.capacity()*sizeof(T);
}

template <typename Key, typename T, class Compare>
long
amrex::bytesOf (const std::map<Key,T,Compare>& m)
{
    return sizeof(m) + m.size()*(sizeof(Key)+sizeof(T)+gcc_map_node_extra_bytes);
}

extern "C" {
    void* amrex_malloc (std::size_t size);
    void amrex_free (void* p);
}

// hash combiner borrowed from Boost
/*
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
template<typename T>
void
amrex::hash_combine (uint64_t & seed, const T & val) noexcept
{
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template<typename T>
uint64_t
amrex::hash_vector (const Vector<T> & vec, uint64_t seed) noexcept
{
    for (const auto & x: vec) {
        hash_combine(seed, x);
    }
    return seed;
}

#endif /*BL_UTILITY_H*/
