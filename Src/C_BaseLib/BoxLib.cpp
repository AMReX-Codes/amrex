#include <unistd.h>
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
#include <ParallelDescriptor.H>
#include <BLProfiler.H>
#include <Utility.H>

#ifndef BL_AMRPROF
#include <ParmParse.H>
#include <DistributionMapping.H>
#include <FArrayBox.H>
#include <FabArray.H>
#include <MultiFab.H>
#endif

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#ifdef BL_USE_F_BASELIB
#include <MemProfiler_f.H>
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <BLBackTrace.H>
#include <MemPool.H>

#if defined(BL_USE_FORTRAN_MPI) || defined(BL_USE_F_INTERFACES)
extern "C" {
    void bl_fortran_mpi_comm_init (int fcomm);
    void bl_fortran_mpi_comm_free ();
    void bl_fortran_sidecar_mpi_comm_free (int fcomm);
}
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

namespace BoxLib
{
    std::string exename;
    int verbose = 0;
}

//
// This is used by BoxLib::Error(), BoxLib::Abort(), and BoxLib::Assert()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//

void
BoxLib::write_to_stderr_without_buffering (const char* str)
{
    //
    // Flush all buffers.
    //
    fflush(NULL);

    if (str)
    {
	std::ostringstream procall;
	procall << ParallelDescriptor::MyProcAll() << "::";
	const char *cprocall = procall.str().c_str();
        const char * const end = " !!!\n";
	fwrite(cprocall, strlen(cprocall), 1, stderr);
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

void
BL_this_is_a_dummy_routine_to_force_version_into_executable ()
{
    BoxLib::write_to_stderr_without_buffering(version);    
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
    BLBackTrace::handler(SIGABRT);
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
    BLBackTrace::handler(SIGABRT);
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

    BLBackTrace::handler(SIGABRT);
}

namespace
{
    std::stack<BoxLib::PTR_TO_VOID_FUNC> The_Finalize_Function_Stack;
    std::stack<BoxLib::PTR_TO_VOID_FUNC> The_Initialize_Function_Stack;
}

void
BoxLib::ExecOnFinalize (PTR_TO_VOID_FUNC fp)
{
    The_Finalize_Function_Stack.push(fp);
}

void
BoxLib::ExecOnInitialize (PTR_TO_VOID_FUNC fp)
{
    The_Initialize_Function_Stack.push(fp);
}

void
BoxLib::Initialize (int& argc, char**& argv, bool build_parm_parse, MPI_Comm mpi_comm)
{
    ParallelDescriptor::StartParallel(&argc, &argv, mpi_comm);

#ifndef WIN32
    //
    // Make sure to catch new failures.
    //
    std::set_new_handler(BoxLib::OutOfMemory);

    if (argv[0][0] != '/') {
	int bufSize(1024);
	char temp[bufSize];
	char *rCheck = getcwd(temp, bufSize);
	if(rCheck == 0) {
	  BoxLib::Abort("**** Error:  getcwd buffer too small.");
	}
	exename = temp;
	exename += "/";
    }
    exename += argv[0];
#endif

#ifdef BL_USE_UPCXX
    upcxx::init(&argc, &argv);
    if (upcxx::myrank() != ParallelDescriptor::MyProc())
	BoxLib::Abort("UPC++ rank != MPI rank");
#endif

#ifdef BL_USE_MPI3
    BL_MPI_REQUIRE( MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &ParallelDescriptor::cp_win) );
    BL_MPI_REQUIRE( MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &ParallelDescriptor::fb_win) );
    BL_MPI_REQUIRE( MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &ParallelDescriptor::fpb_win) );
#endif

    while ( ! The_Initialize_Function_Stack.empty())
    {
        //
        // Call the registered function.
        //
        (*The_Initialize_Function_Stack.top())();
        //
        // And then remove it from the stack.
        //
        The_Initialize_Function_Stack.pop();
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

    signal(SIGSEGV, BLBackTrace::handler); // catch seg falult
    signal(SIGINT,  BLBackTrace::handler);

#ifndef BL_AMRPROF
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

    {
	ParmParse pp("boxlib");
	pp.query("v", verbose);
	pp.query("verbose", verbose);

	int invalid = 0, divbyzero=0, overflow=0;
	pp.query("fpe_trap_invalid", invalid);
	pp.query("fpe_trap_zero", divbyzero);
	pp.query("fpe_trap_overflow", overflow);
	int flags = 0;
	if (invalid)   flags |= FE_INVALID;
	if (divbyzero) flags |= FE_DIVBYZERO;
	if (overflow)  flags |= FE_OVERFLOW;
#if defined(__linux__)
#if !defined(__PGI) || (__PGIC__ >= 16)
	if (flags != 0) {
	    feenableexcept(flags);  // trap floating point exceptions
	    signal(SIGFPE,  BLBackTrace::handler);
	}
#endif
#endif
    }

    ParallelDescriptor::StartTeams();

    ParallelDescriptor::StartSubCommunicator();

    mempool_init();

#endif

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

#if defined(BL_USE_FORTRAN_MPI) || defined(BL_USE_F_INTERFACES)
    int fcomm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
    bl_fortran_mpi_comm_init (fcomm);
#endif

#if defined(BL_MEM_PROFILING) && defined(BL_USE_F_BASELIB)
    MemProfiler_f::initialize();
#endif
}

void
BoxLib::Finalize (bool finalize_parallel)
{
    BL_PROFILE_FINALIZE();

#ifdef BL_LAZY
    Lazy::Finalize();
#endif

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

    // The MemPool stuff is not using The_Finalize_Function_Stack so that
    // it can be used in Fortran BoxLib.
#ifndef BL_AMRPROF
    if (BoxLib::verbose)
    {
	int mp_min, mp_max, mp_tot;
	mempool_get_stats(mp_min, mp_max, mp_tot);  // in MB
	if (ParallelDescriptor::NProcs() == 1) {
	    if (mp_tot > 0) {
		std::cout << "MemPool: " 
#ifdef _OPENMP
			  << "min used in a thread: " << mp_min << " MB, "
			  << "max used in a thread: " << mp_max << " MB, "
#endif
			  << "tot used: " << mp_tot << " MB." << std::endl;
	    }
	} else {
	    int global_max = mp_tot;
	    int global_min = mp_tot;
	    ParallelDescriptor::ReduceIntMax(global_max);
	    if (global_max > 0) {
		ParallelDescriptor::ReduceIntMin(global_min);
		if (ParallelDescriptor::IOProcessor()) {
		    std::cout << "MemPool: " 
			      << "min used in a rank: " << global_min << " MB, "
			      << "max used in a rank: " << global_max << " MB."
			      << std::endl;
		}
	    }
	}
    }
#endif

#ifdef BL_MEM_PROFILING
    MemProfiler::report("Final");
#endif
    
    ParallelDescriptor::EndTeams();

    ParallelDescriptor::EndSubCommunicator();

#ifdef BL_USE_UPCXX
    upcxx::finalize();
#endif

    if (finalize_parallel) {
#if defined(BL_USE_FORTRAN_MPI) || defined(BL_USE_F_INTERFACES)
	bl_fortran_mpi_comm_free();
#endif
    /* Don't shut down MPI if GASNet is still using MPI */
#ifndef GASNET_CONDUIT_MPI
        ParallelDescriptor::EndParallel();
#endif
    }
}

