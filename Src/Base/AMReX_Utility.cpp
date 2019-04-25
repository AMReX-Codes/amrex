#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>
#include <random>

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <AMReX_BLFort.H>
#include <AMReX_REAL.H>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_BLassert.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Print.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/param.h>
#include <unistd.h>


using std::ostringstream;


namespace {
    const char* path_sep_str = "/";
}

//
// Return true if argument is a non-zero length string of digits.
//

bool
amrex::is_integer (const char* str)
{
    int len = 0;

    if (str == 0 || (len = strlen(str)) == 0)
        return false;

    for (int i = 0; i < len; i++)
        if (!isdigit(str[i]))
            return false;

    return true;
}

namespace {
    bool tokenize_initialized = false;
    char* line = 0;
    void CleanupTokenizeStatics ()
    {
        delete [] line;
    }
}

const std::vector<std::string>&
amrex::Tokenize (const std::string& instr,
                  const std::string& separators)
{
    if (!tokenize_initialized) {
        amrex::ExecOnFinalize(CleanupTokenizeStatics);
        tokenize_initialized = true;
    }

    static std::vector<char*>       ptr;
    static std::vector<std::string> tokens;
    static int                      linelen = 0;
    //
    // Make copy of line that we can modify.
    //
    const int len = instr.size() + 1;

    if (len > linelen)
    {
        delete [] line;
        line = new char[len];
        linelen = len;
    }

    (void) std::strcpy(line, instr.c_str());

    char* token = 0;

    if (!((token = std::strtok(line, separators.c_str())) == 0))
    {
        ptr.push_back(token);
        while (!((token = strtok(0, separators.c_str())) == 0))
            ptr.push_back(token);
    }

    tokens.resize(ptr.size());

    for (int i = 1, N = ptr.size(); i < N; i++)
    {
        char* p = ptr[i];

        while (strchr(separators.c_str(), *(p-1)) != 0)
            *--p = 0;
    }

    for (int i = 0, N = ptr.size(); i < N; i++)
        tokens[i] = ptr[i];

    ptr.clear();
    return tokens;
}

std::string
amrex::toLower (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::string
amrex::toUpper (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}


std::string
amrex::Concatenate (const std::string& root,
                     int                num,
                     int                mindigits)
{
    BL_ASSERT(mindigits >= 0);
    std::stringstream result;
    result << root << std::setfill('0') << std::setw(mindigits) << num;
    return result.str();
}


bool
amrex::UtilCreateDirectory (const std::string& path,
			    mode_t mode, bool verbose)
{
    bool retVal(false);
    Vector<std::pair<std::string, int> > pathError;

    if (path.length() == 0 || path == path_sep_str) {
        return true;
    }

    errno = 0;

    if(strchr(path.c_str(), *path_sep_str) == 0) {
        //
        // No slashes in the path.
        //
        errno = 0;
        if(mkdir(path.c_str(), mode) < 0 && errno != EEXIST) {
	  retVal = false;
	} else {
	  retVal = true;
	}
	pathError.push_back(std::make_pair(path, errno));
    } else {
        //
        // Make copy of the directory pathname so we can write to it.
        //
        char *dir = new char[path.length() + 1];
        (void) strcpy(dir, path.c_str());

        char *slash = strchr(dir, *path_sep_str);

        if(dir[0] == *path_sep_str) {  // full pathname.
            do {
                if(*(slash+1) == 0) {
                    break;
		}
                if((slash = strchr(slash+1, *path_sep_str)) != 0) {
                    *slash = 0;
		}
                errno = 0;
                if(mkdir(dir, mode) < 0 && errno != EEXIST) {
	          retVal = false;
		} else {
	          retVal = true;
		}
	        pathError.push_back(std::make_pair(dir, errno));
                if(slash) {
                  *slash = *path_sep_str;
		}
            } while(slash);

        } else {  // relative pathname.

            do {
                *slash = 0;
                errno = 0;
                if(mkdir(dir, mode) < 0 && errno != EEXIST) {
	          retVal = false;
		} else {
	          retVal = true;
		}
	        pathError.push_back(std::make_pair(dir, errno));
                *slash = *path_sep_str;
            } while((slash = strchr(slash+1, *path_sep_str)) != 0);

            errno = 0;
            if(mkdir(dir, mode) < 0 && errno != EEXIST) {
	      retVal = false;
	    } else {
	      retVal = true;
	    }
	    pathError.push_back(std::make_pair(dir, errno));
        }

        delete [] dir;
    }

    if(retVal == false  || verbose == true) {
      for(int i(0); i < pathError.size(); ++i) {
          amrex::AllPrint()<< "amrex::UtilCreateDirectory:: path errno:  " 
                           << pathError[i].first << " :: "
                           << strerror(pathError[i].second)
                           << std::endl;
      }
    }

    return retVal;
}

void
amrex::CreateDirectoryFailed (const std::string& dir)
{
    std::string msg("Couldn't create directory: ");
    msg += dir;
    amrex::Error(msg.c_str());
}

void
amrex::FileOpenFailed (const std::string& file)
{
    std::string msg("Couldn't open file: ");
    msg += file;
    amrex::Error(msg.c_str());
}

void
amrex::UnlinkFile (const std::string& file)
{
    unlink(file.c_str());
}

bool
amrex::FileExists(const std::string &filename)
{
  struct stat statbuff;
  return(::lstat(filename.c_str(), &statbuff) != -1);
}

std::string
amrex::UniqueString()
{
  std::stringstream tempstring;
  tempstring << std::setprecision(11) << std::fixed << ParallelDescriptor::second();
  int tsl(tempstring.str().length());
  return(tempstring.str().substr(tsl/2, tsl));
}

void
amrex::UtilCreateCleanDirectory (const std::string &path, bool callbarrier)
{
  if(ParallelContext::IOProcessorSub()) {
    if(amrex::FileExists(path)) {
      std::string newoldname(path + ".old." + amrex::UniqueString());
      if (amrex::system::verbose > 1) {
          amrex::Print() << "amrex::UtilCreateCleanDirectory():  " << path
                         << " exists.  Renaming to:  " << newoldname << std::endl;
      }
      std::rename(path.c_str(), newoldname.c_str());
    }
    if( ! amrex::UtilCreateDirectory(path, 0755)) {
      amrex::CreateDirectoryFailed(path);
    }
  }
  if(callbarrier) {
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier("amrex::UtilCreateCleanDirectory");
  }
}


void
amrex::UtilCreateDirectoryDestructive(const std::string &path, bool callbarrier)
{
  if(ParallelContext::IOProcessorSub()) 
  {
    if(amrex::FileExists(path)) 
    {
      if (amrex::Verbose() > 1) {
          amrex::Print() << "amrex::UtilCreateCleanDirectoryDestructive():  " << path
                         << " exists.  I am destroying it.  " << std::endl;
      }
      char command[2000];
      sprintf(command, "\\rm -rf %s", path.c_str());;
      int retVal = std::system(command);
      if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
          amrex::Error("Removing old directory failed.");
      }
    }
    if( ! amrex::UtilCreateDirectory(path, 0755)) 
    {
      amrex::CreateDirectoryFailed(path);
    }
  }
  if(callbarrier) 
  {
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier("amrex::UtilCreateCleanDirectoryDestructive");
  }
}

void
amrex::UtilRenameDirectoryToOld (const std::string &path, bool callbarrier)
{
  if(ParallelContext::IOProcessorSub()) {
    if(amrex::FileExists(path)) {
      std::string newoldname(path + ".old." + amrex::UniqueString());
      if (amrex::Verbose() > 1) {
          amrex::Print() << "amrex::UtilRenameDirectoryToOld():  " << path
                         << " exists.  Renaming to:  " << newoldname << std::endl;
      }
      std::rename(path.c_str(), newoldname.c_str());
    }
  }
  if(callbarrier) {
    // Force other processors to wait until directory is renamed.
    ParallelDescriptor::Barrier("amrex::UtilRenameDirectoryToOld");
  }
}

void
amrex::OutOfMemory ()
{
    amrex::Error("Sorry, out of memory, bye ...");
}

namespace
{
    int nthreads;

    amrex::Vector<std::mt19937> generators;

#ifdef AMREX_USE_CUDA
    /**
    * \brief The random seed array is allocated with an extra buffer space to 
    *        reduce the computational cost of dynamic memory allocation and 
    *        random seed generation. 
    */
    __device__ curandState_t *glo_RandStates;
    amrex::Gpu::DeviceVector<curandState_t> dev_RandStates_Seed;
#endif

}

void
amrex::InitRandom (unsigned long seed, int nprocs)
{

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif
    generators.resize(nthreads);

#ifdef _OPENMP
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned long init_seed = seed + tid*nprocs;
        generators[tid].seed(init_seed);
    }
#else
    generators[0].seed(seed);
#endif
}

void amrex::ResetRandomSeed(unsigned long seed)
{
    InitRandom(seed);
}

AMREX_GPU_HOST_DEVICE double
amrex::RandomNormal (double mean, double stddev)
{

    double rand;

#ifdef __CUDA_ARCH__

    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y)) 
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    rand = stddev * curand_normal_double(&glo_RandStates[tid]) + mean; 

#else

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    std::normal_distribution<double> distribution(mean, stddev);
    rand = distribution(generators[tid]);

#endif

    return rand;
}

AMREX_GPU_HOST_DEVICE double
amrex::Random ()
{
    double rand;

#ifdef __CUDA_ARCH__

    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y)) 
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    rand = curand_uniform_double(&glo_RandStates[tid]); 


#else

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    rand = distribution(generators[tid]);

#endif

    return rand;
}


unsigned long
amrex::Random_int(unsigned long n)
{
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    std::uniform_int_distribution<unsigned long> distribution(0, n-1);
    return distribution(generators[tid]);
}

void
amrex::SaveRandomState(std::ostream& os)
{
    for (int i = 0; i < nthreads; i++) {
        os << generators[i] << "\n";
    }
}

void
amrex::RestoreRandomState(std::istream& is, int nthreads_old, int nstep_old)
{
    int N = std::min(nthreads, nthreads_old);
    for (int i = 0; i < N; i++)
        is >> generators[i];
    if (nthreads > nthreads_old) {
        const int NProcs = ParallelDescriptor::NProcs();
        const int MyProc = ParallelDescriptor::MyProc();
        for (int i = nthreads_old; i < nthreads; i++) {
	    unsigned long seed = MyProc+1 + i*NProcs;
	    if (ULONG_MAX/(unsigned long)(nstep_old+1) >static_cast<unsigned long>(nthreads*NProcs)) { // avoid overflow
		seed += nstep_old*nthreads*NProcs;
	    }

            generators[i].seed(seed);
        }
    }
}

void
amrex::UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                           bool printSet)
{
  if(setSize > poolSize) {
    amrex::Abort("**** Error in UniqueRandomSubset:  setSize > poolSize.");
  }
  std::set<int> copySet;
  Vector<int> uSetTemp;
  while(static_cast<int>(copySet.size()) < setSize) {
    int r(amrex::Random_int(poolSize));
    if(copySet.find(r) == copySet.end()) {
      copySet.insert(r);
      uSetTemp.push_back(r);
    }
  }
  uSet = uSetTemp;
  if(printSet) {
    for(int i(0); i < uSet.size(); ++i) {
        amrex::AllPrint() << "uSet[" << i << "]  = " << uSet[i] << std::endl;
    }
  }
}


void 
amrex::InitRandSeedOnDevice (int N)
{
  ResizeRandomSeed(N);
}

void 
amrex::CheckSeedArraySizeAndResize (int N)
{
#ifdef AMREX_USE_CUDA
  if ( dev_RandStates_Seed.size() < N) {
     ResizeRandomSeed(N);
  }
#endif
}

void 
amrex::ResizeRandomSeed (int N)
{

#ifdef AMREX_USE_CUDA  

  int Nbuffer = N * 2;

  int PrevSize = dev_RandStates_Seed.size();

  const int MyProc = amrex::ParallelDescriptor::MyProc();
  int SizeDiff = Nbuffer - PrevSize;

  dev_RandStates_Seed.resize(Nbuffer);
  curandState_t *d_RS_Seed = dev_RandStates_Seed.dataPtr();
  cudaMemcpyToSymbol(glo_RandStates,&d_RS_Seed,sizeof(curandState_t *));

  AMREX_PARALLEL_FOR_1D (SizeDiff, idx,
  {
     unsigned long seed = MyProc*1234567UL + 12345UL ;
     int seqstart = idx + 10 * idx ; 
     int loc = idx + PrevSize;
     curand_init(seed, seqstart, 0, &glo_RandStates[loc]);
  }); 

#endif

}

void 
amrex::DeallocateRandomSeedDevArray()
{
#ifdef AMREX_USE_CUDA  
  dev_RandStates_Seed.resize(0);
  dev_RandStates_Seed.shrink_to_fit();
#endif
}



void
amrex::NItemsPerBin (int totalItems, Vector<int> &binCounts)
{
  if(binCounts.size() == 0) {
    return;
  }
  bool verbose(false);
  int countForAll(totalItems / binCounts.size());
  int remainder(totalItems % binCounts.size());
  if(verbose) {
      amrex::Print() << "amrex::NItemsPerBin:  countForAll remainder = " << countForAll
                     << "  " << remainder << std::endl;
  }
  for(int i(0); i < binCounts.size(); ++i) {
    binCounts[i] = countForAll;
  }
  for(int i(0); i < remainder; ++i) {
    ++binCounts[i];
  }
  for(int i(0); i < binCounts.size(); ++i) {
    if(verbose) {
        amrex::Print() << "amrex::NItemsPerBin::  binCounts[" << i << "] = " << binCounts[i] << std::endl;
    }
  }
}

// -------------------------------------------------------------------
int amrex::CRRBetweenLevels(int fromlevel, int tolevel,
                            const Vector<int> &refratios)
{
  BL_ASSERT(fromlevel >= 0);
  BL_ASSERT(tolevel >= fromlevel);
  BL_ASSERT(tolevel <= refratios.size());
  int level, rr = 1;
  for(level = fromlevel; level < tolevel; ++level) {
    rr *= refratios[level];
  }
  return rr;
}

//
// Fortran entry points for amrex::Random().
//

#ifndef AMREX_XSDK
BL_FORT_PROC_DECL(BLUTILINITRAND,blutilinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    amrex::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLINITRAND,blinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    amrex::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLUTILRAND,blutilrand)(amrex::Real* rn)
{
    *rn = amrex::Random();
}
#endif

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
amrex::InvNormDist (double p)
{
    if (p <= 0 || p >= 1)
        amrex::Error("amrex::InvNormDist(): p MUST be in (0,1)");
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

#ifndef AMREX_XSDK
BL_FORT_PROC_DECL(BLINVNORMDIST,blinvnormdist)(amrex::Real* result)
{
    //
    // Convert from [0, 1) to (0,1)
    // 
    double val = 0.0;
    while (val == 0.0) {
        val = amrex::Random();
    }

    *result = amrex::InvNormDist(val);
}
#endif


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
amrex::InvNormDistBest (double p)

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
      amrex::Error("InvNormDistBest(): p MUST be in (0,1)");

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

#ifndef AMREX_XSDK
BL_FORT_PROC_DECL(BLINVNORMDISTBEST,blinvnormdistbest)(amrex::Real* result)
{
    //
    // Convert from [0, 1) to (0,1)
    // 
    double val = 0.0;
    while (val == 0.0) {
        val = amrex::Random();
    }

    *result = amrex::InvNormDistBest(val);
}
#endif

//
// Sugar for parsing IO
//

std::istream&
amrex::operator>>(std::istream& is, const expect& exp)
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
	amrex::Error(msg.c_str());
    }
    return is;
}

amrex::expect::expect(const char* istr_)
    : istr(istr_)
{
}

amrex::expect::expect(const std::string& str_)
    : istr(str_)
{
}

amrex::expect::expect(char c)
{
    istr += c;
}

const std::string&
amrex::expect::the_string() const
{
    return istr;
}


//
// StreamRetry
//

int amrex::StreamRetry::nStreamErrors = 0;

amrex::StreamRetry::StreamRetry(std::ostream &a_os, const std::string &a_suffix,
                                 const int a_maxtries)
    : tries(0), maxTries(a_maxtries), sros(a_os), spos(a_os.tellp()), suffix(a_suffix)
{
}

amrex::StreamRetry::StreamRetry(const std::string &filename,
				 const bool abortonretryfailure,
                                 const int maxtries)
    : tries(0), maxTries(maxtries),
      abortOnRetryFailure(abortonretryfailure),
      fileName(filename),
      sros(amrex::ErrorStream())    // unused here, just to make the compiler happy
{
  nStreamErrors = 0;
}

bool amrex::StreamRetry::TryOutput()
{
  if(tries == 0) {
    ++tries;
    return true;
  } else {
    if(sros.fail()) {
      ++nStreamErrors;
      int myProc(ParallelDescriptor::MyProc());
      if(tries <= maxTries) {
          if (amrex::Verbose() > 1) {
              amrex::AllPrint() << "PROC: " << myProc << " :: STREAMRETRY_" << suffix << " # "
                                << tries << " :: gbfe:  "
                                << sros.good() << sros.bad() << sros.fail() << sros.eof()
                                << " :: sec = " << ParallelDescriptor::second()
                                << " :: os.tellp() = " << sros.tellp()
                                << " :: rewind spos = " << spos
                                << std::endl;
          }
        sros.clear();  // clear the bad bits
        if (amrex::Verbose() > 1) {
            amrex::AllPrint() << "After os.clear() : gbfe:  " << sros.good() << sros.bad()
                              << sros.fail() << sros.eof() << std::endl;
        }
        sros.seekp(spos, std::ios::beg);  // reset stream position
        ++tries;
        return true;
      } else {
        if (amrex::Verbose() > 1) {
            amrex::AllPrint() << "PROC: " << myProc << " :: STREAMFAILED_" << suffix << " # "
                              << tries << " :: File may be corrupt.  :: gbfe:  "
                              << sros.good() << sros.bad() << sros.fail() << sros.eof()
                              << " :: sec = " << ParallelDescriptor::second()
                              << " :: os.tellp() = " << sros.tellp()
                              << " :: rewind spos = " << spos
                              << std::endl;
        }
        sros.clear();  // clear the bad bits
        if (amrex::Verbose() > 1) {
            amrex::AllPrint() << "After os.clear() : gbfe:  " << sros.good() << sros.bad()
                              << sros.fail() << sros.eof() << std::endl;
        }
        return false;
      }
    } else {
      return false;
    }
  }
}


bool amrex::StreamRetry::TryFileOutput()
{
    bool bTryOutput(false);

    if(tries == 0) {
      bTryOutput = true;
    } else {

      int nWriteErrors(nStreamErrors);
      ParallelDescriptor::ReduceIntSum(nWriteErrors);

      if(nWriteErrors == 0) {  // wrote a good file
        bTryOutput = false;
      } else {                 // wrote a bad file, rename it
        if(ParallelDescriptor::IOProcessor()) {
          const std::string& badFileName = amrex::Concatenate(fileName + ".bad",
                                                               tries - 1, 2);
          if (amrex::Verbose() > 1) {
              amrex::Print() << nWriteErrors << " STREAMERRORS : Renaming file from "
                             << fileName << "  to  " << badFileName << std::endl;
          }
          std::rename(fileName.c_str(), badFileName.c_str());
        }
        ParallelDescriptor::Barrier("StreamRetry::TryFileOutput");  // wait for file rename

        // check for maxtries and abort pref
        if(tries < maxTries) {
          bTryOutput = true;
        } else {
          if(abortOnRetryFailure) {
            amrex::Abort("STREAMERROR : StreamRetry::maxTries exceeded.");
          }
          bTryOutput = false;
        }
      }
    }

    ++tries;
    nStreamErrors = 0;
    return bTryOutput;
}


void amrex::SyncStrings(const Vector<std::string> &localStrings,
                         Vector<std::string> &syncedStrings, bool &alreadySynced)
{
#ifdef BL_USE_MPI
  const int nProcs(ParallelDescriptor::NProcs());
  const int ioProcNumber(ParallelDescriptor::IOProcessorNumber());
  int nUnmatched(0);

  Vector<std::string> localStringsCopy = localStrings;

  // ---- broadcast ioproc strings
  int pfStringsSize(0);
  std::ostringstream pfStrings;
  if(ParallelDescriptor::IOProcessor()) {
    for(int i(0); i < localStringsCopy.size(); ++i) {
      pfStrings << localStringsCopy[i] << '\n';
    }
    pfStringsSize = pfStrings.str().size();
  }
  ParallelDescriptor::Bcast(&pfStringsSize, 1);

  Vector<char> pfCharArray(pfStringsSize + 1);
  if(ParallelDescriptor::IOProcessor()) {
    std::strcpy(pfCharArray.dataPtr(), pfStrings.str().c_str());  // null terminated
  }
  ParallelDescriptor::Bcast(pfCharArray.dataPtr(), pfCharArray.size());

  // ---- extract the ioproc strings
  Vector<std::string> ioprocStrings, sendStrings;
  if( ! ParallelDescriptor::IOProcessor()) {
    std::istringstream pfIn(pfCharArray.dataPtr());
    std::string pfName;
    while( ! pfIn.eof()) {
      std::getline(pfIn, pfName, '\n');
      if( ! pfIn.eof()) {
	ioprocStrings.push_back(pfName);
      }
    }
    // ---- now check if they match on non ioprocs
    for(int n(0); n < ioprocStrings.size(); ++n) {
      bool matched(false);
      for(int i(0); i < localStringsCopy.size(); ++i) {
        if(ioprocStrings[n] == localStringsCopy[i]) {
          matched = true;
        }
      }
      if( ! matched) {
        ++nUnmatched;
        localStringsCopy.push_back(ioprocStrings[n]);  // ---- add to local set
      }
    }
    for(int n(0); n < localStringsCopy.size(); ++n) {
      bool matched(false);
      for(int i(0); i < ioprocStrings.size(); ++i) {
        if(localStringsCopy[n] == ioprocStrings[i]) {
	  matched = true;
        }
      }
      if( ! matched) {
        ++nUnmatched;
        sendStrings.push_back(localStringsCopy[n]);  // ---- send these to the ioproc
      }
    }
  }

  ParallelDescriptor::ReduceIntMax(nUnmatched);
  if(nUnmatched == 0) {
    alreadySynced = true;
    syncedStrings = localStrings;
    return;
  }
  alreadySynced = false;


  int sendStringsSize(0);
  std::ostringstream ossSendStrings;
  Vector<char> sendCharArray(1);  // cannot be zero for gather call
  if( ! ParallelDescriptor::IOProcessor()) {
    for(int i(0); i < sendStrings.size(); ++i) {
      ossSendStrings << sendStrings[i] << '\n';
    }
    sendStringsSize = ossSendStrings.str().size();
    sendCharArray.resize(sendStringsSize + 1);
    std::strcpy(sendCharArray.dataPtr(), ossSendStrings.str().c_str());  // null terminated
  }

  Vector<int> nChars(nProcs, 0);
  ParallelDescriptor::Gather(&sendStringsSize, 1, nChars.dataPtr(), 1, ioProcNumber);

  int totalChars(0);
  Vector<char> recvStrings(1);
  Vector<int> offset(nProcs, 0);
  if(ParallelDescriptor::IOProcessor()) {
    for(int i(0); i < nChars.size(); ++i) {
      totalChars += nChars[i];
    }
    recvStrings.resize(totalChars + 1);

    int iOffset(0);
    for(int i(1); i < nChars.size(); ++i) {
      iOffset += nChars[i-1];
      offset[i] = iOffset;
    }
  }

  BL_MPI_REQUIRE( MPI_Gatherv(sendCharArray.dataPtr(),
                              sendStringsSize, 
                              ParallelDescriptor::Mpi_typemap<char>::type(),
                              recvStrings.dataPtr(),
                              nChars.dataPtr(),
                              offset.dataPtr(),
                              ParallelDescriptor::Mpi_typemap<char>::type(),
                              ioProcNumber,
                              ParallelDescriptor::Communicator()) );

  if(ParallelDescriptor::IOProcessor()) {
    std::istringstream pfIn(recvStrings.dataPtr());
    std::string pfName;
    syncedStrings = localStrings;
    while( ! pfIn.eof()) {
      std::getline(pfIn, pfName, '\n');
      if( ! pfIn.eof()) {
	syncedStrings.push_back(pfName);  // ---- add the gathered strings
      }
    }
  }

  // ---- broadcast synced ioproc strings
  int syncedStringsSize(0);
  std::ostringstream syncedStrStr;
  if(ParallelDescriptor::IOProcessor()) {
    for(int i(0); i < syncedStrings.size(); ++i) {
      syncedStrStr << syncedStrings[i] << '\n';
    }
    syncedStringsSize = syncedStrStr.str().size();
  }
  ParallelDescriptor::Bcast(&syncedStringsSize, 1);

  Vector<char> syncedCharArray(syncedStringsSize + 1);
  if(ParallelDescriptor::IOProcessor()) {
    std::strcpy(syncedCharArray.dataPtr(), syncedStrStr.str().c_str());  // null terminated
  }
  ParallelDescriptor::Bcast(syncedCharArray.dataPtr(), syncedCharArray.size());

  if( ! ParallelDescriptor::IOProcessor()) {
    std::istringstream syncedIn(syncedCharArray.dataPtr());
    std::string syncedName;
    while( ! syncedIn.eof()) {
      std::getline(syncedIn, syncedName, '\n');
      if( ! syncedIn.eof()) {
	syncedStrings.push_back(syncedName);
      }
    }
  }

#else
    alreadySynced = true;
    syncedStrings = localStrings;
#endif

}

amrex::Vector<char> amrex::SerializeStringArray(const Vector<std::string> &stringArray)
{
  std::ostringstream stringStream;
  for(int i(0); i < stringArray.size(); ++i) {
    stringStream << stringArray[i] << '\n';
  }

  Vector<char> charArray(stringStream.str().size() + 1);
  std::strcpy(charArray.dataPtr(), stringStream.str().c_str());  // null terminated

  return charArray;
}

amrex::Vector<std::string> amrex::UnSerializeStringArray(const Vector<char> &charArray)
{
  Vector<std::string> stringArray;
  std::istringstream stringStream(charArray.dataPtr());
  std::string sTemp;
  while( ! stringStream.eof()) {
    std::getline(stringStream, sTemp, '\n');
    if( ! stringStream.eof()) {
      stringArray.push_back(sTemp);
    }
  }

  return stringArray;
}

void amrex::BroadcastBool(bool &bBool, int myLocalId, int rootId, const MPI_Comm &localComm)
{
  int numBool;
  if (myLocalId == rootId) {
    numBool = bBool;
  }

  amrex::ParallelDescriptor::Bcast(&numBool, 1, rootId, localComm);

  if(myLocalId != rootId) {
    bBool = numBool;
  }
}


void amrex::BroadcastString(std::string &bStr, int myLocalId, int rootId, const MPI_Comm &localComm)
{
  Vector<std::string> vecString(1, bStr);
  Vector<char> serialString;
  if(myLocalId == rootId) {
    serialString = amrex::SerializeStringArray(vecString);
  }

  amrex::BroadcastArray(serialString, myLocalId, rootId, localComm);

  if(myLocalId != rootId) {
    vecString = amrex::UnSerializeStringArray(serialString);
    bStr = vecString[0];
  }
}

void amrex::BroadcastStringArray(Vector<std::string> &bSA, int myLocalId, int rootId, const MPI_Comm &localComm)
{
  Vector<char> serialStringArray;
  if(myLocalId == rootId) {
    serialStringArray = amrex::SerializeStringArray(bSA);
  }

  amrex::BroadcastArray(serialStringArray, myLocalId, rootId, localComm);

  if(myLocalId != rootId) {
    bSA = amrex::UnSerializeStringArray(serialStringArray);
  }
}

void amrex::USleep(double sleepsec) {
  constexpr unsigned int msps = 1000000;
  usleep(sleepsec * msps);
}


namespace {
    static auto clock_time_begin = amrex::MaxResSteadyClock::now();
}

double amrex::second () noexcept
{
    return std::chrono::duration_cast<std::chrono::duration<double> >
        (amrex::MaxResSteadyClock::now() - clock_time_begin).count();
}


extern "C" {
    void* amrex_malloc (std::size_t size)
    {
        return malloc(size);
    }


    void amrex_free (void* p)
    {
        std::free(p);
    }

    double amrex_random ()
    {
        return amrex::Random();
    }

    long amrex_random_int (long n)  // This is for Fortran, which doesn't have unsigned long.
    {
        return static_cast<long>(amrex::Random_int(static_cast<unsigned long>(n)));
    }
}
