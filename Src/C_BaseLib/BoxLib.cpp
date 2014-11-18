
#include <winstd.H>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <new>
#include <stack>
#include <limits>

#include <BoxLib.H>
#include <DistributionMapping.H>
#include <FArrayBox.H>
#include <FabArray.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <Profiler.H>
#include <Utility.H>
#include <MultiFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#define bl_str(s)  # s
#define bl_xstr(s) bl_str(s)
//
// The definition of our version string.
//    
// Takes the form:  boxlib version 2.0 built Jun 25 1996 at 14:52:36
//
const char * const version =

"boxlib version "
bl_xstr(BL_VERSION_MAJOR)
"."
bl_xstr(BL_VERSION_MINOR)
" built "
__DATE__
" at "
__TIME__;

#undef bl_str
#undef bl_xstr

//
// This is used by BoxLib::Error(), BoxLib::Abort(), and BoxLib::Assert()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//

static
void
write_to_stderr_without_buffering (const char* str)
{
    //
    // Flush all buffers.
    //
    fflush(NULL);

    if (str)
    {
        const char * const end = " !!!\n";
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

void
BL_this_is_a_dummy_routine_to_force_version_into_executable ()
{
    write_to_stderr_without_buffering(version);    
}

static
void
write_lib_id(const char* msg)
{
    fflush(0);
    const char* const boxlib = "BoxLib::";
    fwrite(boxlib, strlen(boxlib), 1, stderr);
    if ( msg ) 
    {
	fwrite(msg, strlen(msg), 1, stderr);
	fwrite("::", 2, 1, stderr);
    }
}

void
BoxLib::Error (const char* msg)
{
    write_lib_id("Error");
    write_to_stderr_without_buffering(msg);
    ParallelDescriptor::Abort();
}

namespace
{
  const int EOS = -1;

  std::string
  Trim (const std::string& str)
  {
    int n;
    for ( n = str.size(); --n >= 0; )
      {
	if ( str[n] != ' ' ) break;
      }
    std::string result;
    for (int i = 0; i <= n; ++i )
      {
	result += str[i];
      }
    return result;
  }

  std::string
  Fint_2_string (const int* iarr, int nlen)
  {
    std::string res;
    for ( int i = 0; i < nlen && *iarr != EOS; ++i )
      {
	res += *iarr++;
      }
    return Trim(res);
  }
}

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_ERROR_CPP,bl_error_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = "FORTRAN:";
  res += Fint_2_string(istr, *NSTR);
  BoxLib::Error(res.c_str());
}

BL_FORT_PROC_DECL(BL_WARNING_CPP,bl_warning_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = "FORTRAN:";
  res += Fint_2_string(istr, *NSTR);
  BoxLib::Warning(res.c_str());
}

BL_FORT_PROC_DECL(BL_ABORT_CPP,bl_abort_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = "FORTRAN:";
  res += Fint_2_string(istr, *NSTR);
  BoxLib::Abort(res.c_str());
}

void
BoxLib::Abort (const char* msg)
{
    write_lib_id("Abort");
    write_to_stderr_without_buffering(msg);
    ParallelDescriptor::Abort();
}

void
BoxLib::Warning (const char* msg)
{
    if (msg)
    {
        std::cerr << msg << '!' << '\n';
    }
}

void
BoxLib::Assert (const char* EX,
                const char* file,
                int         line)
{
    const int N = 512;

    char buf[N];

    snprintf(buf,
             N,
             "Assertion `%s' failed, file \"%s\", line %d",
             EX,
             file,
             line);

    write_to_stderr_without_buffering(buf);

    ParallelDescriptor::Abort();
}

namespace
{
    std::stack<BoxLib::PTR_TO_VOID_FUNC> The_Finalize_Function_Stack;
}

void
BoxLib::ExecOnFinalize (PTR_TO_VOID_FUNC fp)
{
    The_Finalize_Function_Stack.push(fp);
}

void
BoxLib::Initialize (int& argc, char**& argv, bool build_parm_parse, MPI_Comm mpi_comm)
{
#ifndef WIN32
    //
    // Make sure to catch new failures.
    //
    std::set_new_handler(BoxLib::OutOfMemory);
#endif

    ParallelDescriptor::StartParallel(&argc, &argv, mpi_comm);

    if(ParallelDescriptor::NProcsPerfMon() > 0) {
      if(ParallelDescriptor::MyProcAll() == ParallelDescriptor::MyProcAllPerfMon()) {
        std::cout << "Starting PerfMonProc:  myprocall = "
                  << ParallelDescriptor::MyProcAll() << std::endl;
        Profiler::PerfMonProcess();
        BoxLib::Finalize();
        return;
      }
    }

    BL_PROFILE_INITIALIZE();

    //
    // Initialize random seed after we're running in parallel.
    //
    BoxLib::InitRandom(ParallelDescriptor::MyProc()+1, ParallelDescriptor::NProcs());

#ifdef BL_USE_MPI
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "MPI initialized with "
                  << ParallelDescriptor::NProcs()
                  << " MPI processes\n";
    }
#endif

#ifdef _OPENMP
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "OMP initialized with "
                  << omp_get_max_threads()
                  << " OMP threads\n";
    }
#endif

    if (build_parm_parse)
    {
        if (argc == 1)
        {
            ParmParse::Initialize(0,0,0);
        }
        else
        {
            if (strchr(argv[1],'='))
            {
                ParmParse::Initialize(argc-1,argv+1,0);
            }
            else
            {
                ParmParse::Initialize(argc-2,argv+2,argv[1]);
            }
        }
    }

    std::cout << std::setprecision(10);

    if (double(std::numeric_limits<long>::max()) < 9.e18)
    {
	if (ParallelDescriptor::IOProcessor())
	{
	    std::cout << "!\n! WARNING: Maximum of long int, "
		      << std::numeric_limits<long>::max() 
		      << ", might be too small for big runs.\n!\n";
	}
    }
}

void
BoxLib::Finalize (bool finalize_parallel)
{
    BL_PROFILE_FINALIZE();

    while (!The_Finalize_Function_Stack.empty())
    {
        //
        // Call the registered function.
        //
        (*The_Finalize_Function_Stack.top())();
        //
        // And then remove it from the stack.
        //
        The_Finalize_Function_Stack.pop();
    }

    if (finalize_parallel)
        ParallelDescriptor::EndParallel();
}

